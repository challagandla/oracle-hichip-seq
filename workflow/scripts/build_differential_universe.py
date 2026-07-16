"""Build an exact-grid, q-unthresholded differential hypothesis universe.

Only a predeclared raw-count abundance filter is used.  FitHiChIP p/q values and
case/control labels never enter candidate selection.  A temporary SQLite table
keeps multi-million-row inputs bounded in memory and enforces one observation per
sample and exact contact pixel.
"""
import hashlib
import json
import logging
import os
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from normalize_fithichip_all import (  # noqa: E402
    OUT_COLUMNS,
    OUT_COORDS,
    PRIMARY,
    _overlap_mask,
    _read_chromsizes,
    _read_interval_index,
)
from utils import setup_logging  # noqa: E402


log = logging.getLogger(__name__)
SQL_COORDS = ", ".join(OUT_COORDS)
SUPPORT_COLUMNS = OUT_COORDS + [
    "sample_support", "support_samples", "source_counts", "tolerance_bins",
    "source_kind", "min_count", "min_samples",
]


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_audits(
    paths: list[str], expected_samples: list[str], min_count: int,
    bin_size: int, lower_distance: int, upper_distance: int,
) -> dict[str, dict]:
    audits: dict[str, dict] = {}
    for path in paths:
        audit = json.loads(Path(path).read_text(encoding="utf-8"))
        sample = str(audit.get("sample", ""))
        if not sample or sample in audits:
            raise ValueError(f"Duplicate or empty sample in normalization audit: {path}")
        expected = {
            "schema": "oracle-fithichip-all-interactions-v1",
            "source_kind": "fithichip_all_interactions",
            "fithichip_q_filter": None,
            "merge_nearby": False,
            "eligible_min_count": min_count,
            "bin_size": bin_size,
            "lower_distance": lower_distance,
            "upper_distance": upper_distance,
        }
        for key, value in expected.items():
            if audit.get(key) != value:
                raise ValueError(
                    f"Normalization audit {sample!r} has {key}={audit.get(key)!r}; "
                    f"expected {value!r}"
                )
        audits[sample] = audit
    if set(audits) != set(expected_samples):
        raise ValueError(
            "Normalization audit samples must equal comparison samples; "
            f"expected={sorted(expected_samples)}, observed={sorted(audits)}"
        )
    return audits


def _validate_chunk(
    chunk: pd.DataFrame,
    *,
    sample: str,
    bin_size: int,
    lower_distance: int,
    upper_distance: int,
    min_count: int,
    require_min_count: bool,
    chromsizes: dict[str, int],
    blacklist_index: dict[str, tuple[np.ndarray, np.ndarray]],
) -> pd.DataFrame:
    missing = [column for column in OUT_COLUMNS if column not in chunk.columns]
    if missing:
        raise ValueError(f"All-interaction table for {sample} lacks columns: {missing}")
    out = chunk[OUT_COLUMNS].copy()
    for column in ["start1", "end1", "start2", "end2", "score"]:
        numeric = pd.to_numeric(out[column], errors="coerce")
        if numeric.isna().any() or not np.equal(numeric, np.floor(numeric)).all():
            raise ValueError(f"All-interaction table for {sample} has invalid {column}")
        out[column] = numeric.astype(np.int64)
    for column in ["pvalue", "fdr"]:
        numeric = pd.to_numeric(out[column], errors="coerce")
        malformed = out[column].notna() & numeric.isna()
        invalid = numeric.notna() & (
            (~np.isfinite(numeric)) | (numeric < 0) | (numeric > 1)
        )
        if malformed.any() or invalid.any():
            raise ValueError(
                f"All-interaction table for {sample} has invalid {column} values"
            )
        out[column] = numeric.astype(float)
    primary = out["chrom1"].astype(str).str.fullmatch(PRIMARY)
    primary &= out["chrom2"].astype(str).str.fullmatch(PRIMARY)
    distance = out["start2"] - out["start1"]
    count_valid = out["score"] >= min_count if require_min_count else out["score"] > 0
    valid = (
        primary
        & (out["chrom1"].astype(str) == out["chrom2"].astype(str))
        & (out["start1"] >= 0)
        & (out["start2"] >= 0)
        & (out["end1"] - out["start1"] == bin_size)
        & (out["end2"] - out["start2"] == bin_size)
        & (out["start1"] % bin_size == 0)
        & (out["start2"] % bin_size == 0)
        & (out["start1"] < out["start2"])
        & (distance >= lower_distance)
        & (distance <= upper_distance)
        & count_valid
    )
    unknown = ~out["chrom1"].astype(str).isin(chromsizes)
    if unknown.any():
        valid &= ~unknown
    known = ~unknown
    if known.any():
        chrom_end = out.loc[known, "chrom1"].map(chromsizes).astype(np.int64)
        valid.loc[known] &= (
            (out.loc[known, "end1"] <= chrom_end)
            & (out.loc[known, "end2"] <= chrom_end)
        )
    blocked = _overlap_mask(
        out["chrom1"], out["start1"], out["end1"], blacklist_index
    ) | _overlap_mask(
        out["chrom2"], out["start2"], out["end2"], blacklist_index
    )
    valid &= ~blocked
    if not valid.all():
        raise ValueError(
            f"All-interaction table for {sample} has {int((~valid).sum())} rows outside "
            "the exact-grid/count/distance/blacklist contract"
        )
    return out


def _create_database(path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.execute("PRAGMA temp_store=FILE")
    conn.execute(
        """
        CREATE TABLE observations (
            chrom1 TEXT NOT NULL, start1 INTEGER NOT NULL, end1 INTEGER NOT NULL,
            chrom2 TEXT NOT NULL, start2 INTEGER NOT NULL, end2 INTEGER NOT NULL,
            sample TEXT NOT NULL, score INTEGER NOT NULL,
            PRIMARY KEY (chrom1,start1,end1,chrom2,start2,end2,sample)
        ) WITHOUT ROWID
        """
    )
    return conn


def build_universe(
    eligible_files: list[str],
    all_interaction_files: list[str],
    audit_files: list[str],
    blacklist: str | Path,
    chromsizes_path: str | Path,
    bedpe_output: str | Path,
    support_output: str | Path,
    manifest_output: str | Path,
    *,
    comparison: str,
    expected_samples: list[str],
    bin_size: int,
    lower_distance: int,
    upper_distance: int,
    min_count: int,
    min_samples: int,
    chunk_size: int = 250_000,
) -> dict:
    if len(eligible_files) != len(expected_samples):
        raise ValueError("Eligible table count must equal expected sample count")
    if len(all_interaction_files) != len(expected_samples):
        raise ValueError("All-interaction table count must equal expected sample count")
    if len(expected_samples) != len(set(expected_samples)):
        raise ValueError("Comparison sample IDs must be unique")
    if min_count < 1 or min_samples < 1 or min_samples > len(expected_samples):
        raise ValueError("Invalid differential min_count/min_samples contract")
    audits = _load_audits(
        audit_files, expected_samples, min_count, bin_size,
        lower_distance, upper_distance,
    )
    chromsizes = _read_chromsizes(chromsizes_path)
    blacklist_index = _read_interval_index(blacklist)

    manifest_target = Path(manifest_output)
    manifest_target.parent.mkdir(parents=True, exist_ok=True)
    database = manifest_target.with_name(f".{manifest_target.name}.sqlite-{os.getpid()}")
    bed_tmp = Path(bedpe_output).with_name(f".{Path(bedpe_output).name}.tmp-{os.getpid()}")
    support_tmp = Path(support_output).with_name(
        f".{Path(support_output).name}.tmp-{os.getpid()}"
    )
    eligible_rows: dict[str, int] = {}
    all_rows: dict[str, int] = {}
    conn: sqlite3.Connection | None = None
    try:
        conn = _create_database(database)
        insert = (
            "INSERT INTO observations VALUES (?,?,?,?,?,?,?,?)"
        )
        # Pass 1 is deliberately limited to the compact abundance-eligible streams.
        # Only these rows are needed to derive the label-blind candidate key set.
        for path, expected_sample in zip(eligible_files, expected_samples):
            sample = str(audits[expected_sample]["sample"])
            source_path = Path(path)
            if source_path.parent.name != expected_sample or not source_path.name.startswith(
                f"{expected_sample}."
            ):
                raise ValueError(
                    f"Eligible table {path} does not match expected sample "
                    f"{expected_sample!r}"
                )
            observed_rows = 0
            reader = pd.read_csv(path, sep="\t", chunksize=chunk_size)
            try:
                for chunk in reader:
                    out = _validate_chunk(
                        chunk,
                        sample=sample,
                        bin_size=bin_size,
                        lower_distance=lower_distance,
                        upper_distance=upper_distance,
                        min_count=min_count,
                        require_min_count=True,
                        chromsizes=chromsizes,
                        blacklist_index=blacklist_index,
                    )
                    records = [
                        (*coords, sample, int(score))
                        for *coords, score in out[OUT_COORDS + ["score"]].itertuples(
                            index=False, name=None
                        )
                    ]
                    try:
                        conn.executemany(insert, records)
                    except sqlite3.IntegrityError as exc:
                        raise ValueError(
                            f"Eligible table for {sample} contains duplicate contact pixels"
                        ) from exc
                    observed_rows += len(out)
            except pd.errors.EmptyDataError as exc:
                raise ValueError(f"Eligible table for {sample} has no header") from exc
            eligible_rows[sample] = observed_rows
            expected_rows = int(audits[sample]["retained_abundance_eligible_rows"])
            if observed_rows != expected_rows:
                raise ValueError(
                    f"Eligible row count for {sample} is {observed_rows}; audit says "
                    f"{expected_rows}"
                )
            conn.commit()

        conn.execute(
            f"""
            CREATE TEMP TABLE candidates AS
            SELECT {SQL_COORDS}, COUNT(*) AS sample_support
            FROM observations
            GROUP BY {SQL_COORDS}
            HAVING COUNT(*) >= ?
            """,
            (min_samples,),
        )
        candidate_count = int(
            conn.execute("SELECT COUNT(*) FROM candidates").fetchone()[0]
        )
        if candidate_count == 0:
            raise RuntimeError(
                f"No exact-grid contact has cc >= {min_count} in at least "
                f"{min_samples} selected samples"
            )

        candidate_support_by_key = {
            tuple(row[:6]): int(row[6])
            for row in conn.execute(
                f"SELECT {SQL_COORDS}, sample_support FROM candidates"
            )
        }
        candidate_keys = set(candidate_support_by_key)
        conn.execute(
            """
            CREATE TABLE candidate_counts (
                chrom1 TEXT NOT NULL, start1 INTEGER NOT NULL, end1 INTEGER NOT NULL,
                chrom2 TEXT NOT NULL, start2 INTEGER NOT NULL, end2 INTEGER NOT NULL,
                sample TEXT NOT NULL, score INTEGER NOT NULL,
                PRIMARY KEY (chrom1,start1,end1,chrom2,start2,end2,sample)
            ) WITHOUT ROWID
            """
        )
        insert_candidate = "INSERT INTO candidate_counts VALUES (?,?,?,?,?,?,?,?)"

        # Pass 2 streams every q-unthresholded normalized row but retains only exact
        # candidate keys. This captures 1..(min_count-1) source observations while
        # keeping RAM and SQLite proportional to candidates, not 6-10M raw rows/sample.
        for path, expected_sample in zip(all_interaction_files, expected_samples):
            sample = str(audits[expected_sample]["sample"])
            source_path = Path(path)
            if source_path.parent.name != expected_sample or not source_path.name.startswith(
                f"{expected_sample}."
            ):
                raise ValueError(
                    f"All-interaction table {path} does not match expected sample "
                    f"{expected_sample!r}"
                )
            observed_rows = 0
            reader = pd.read_csv(path, sep="\t", chunksize=chunk_size)
            try:
                for chunk in reader:
                    out = _validate_chunk(
                        chunk,
                        sample=sample,
                        bin_size=bin_size,
                        lower_distance=lower_distance,
                        upper_distance=upper_distance,
                        min_count=min_count,
                        require_min_count=False,
                        chromsizes=chromsizes,
                        blacklist_index=blacklist_index,
                    )
                    observed_rows += len(out)
                    records = []
                    for row in out[OUT_COORDS + ["score"]].itertuples(
                        index=False, name=None
                    ):
                        coords, score = tuple(row[:6]), int(row[6])
                        if coords in candidate_keys:
                            records.append((*coords, sample, score))
                    try:
                        conn.executemany(insert_candidate, records)
                    except sqlite3.IntegrityError as exc:
                        raise ValueError(
                            f"All-interaction table for {sample} contains a duplicate "
                            "candidate pixel"
                        ) from exc
            except pd.errors.EmptyDataError as exc:
                raise ValueError(f"All-interaction table for {sample} has no header") from exc
            all_rows[sample] = observed_rows
            expected_rows = int(audits[sample]["retained_all_rows"])
            if observed_rows != expected_rows:
                raise ValueError(
                    f"All-interaction row count for {sample} is {observed_rows}; "
                    f"audit says {expected_rows}"
                )
            conn.commit()

        all_match = " AND ".join(
            [f"a.{column} = e.{column}" for column in OUT_COORDS]
            + ["a.sample = e.sample"]
        )
        inconsistent_eligible = int(conn.execute(
            f"""
            SELECT COUNT(*)
            FROM observations e
            JOIN candidates c USING ({SQL_COORDS})
            LEFT JOIN candidate_counts a ON {all_match}
            WHERE a.score IS NULL OR a.score != e.score
            """
        ).fetchone()[0])
        if inconsistent_eligible:
            raise ValueError(
                f"{inconsistent_eligible} abundance-eligible candidate observations "
                "disagree with their normalized all-interaction source rows"
            )

        query = f"""
            SELECT {', '.join('c.' + c for c in OUT_COORDS)}, o.sample, o.score
            FROM candidates c
            LEFT JOIN candidate_counts o USING ({SQL_COORDS})
            ORDER BY c.chrom1, c.start1, c.end1, c.chrom2, c.start2, c.end2,
                     o.sample
        """
        Path(bedpe_output).parent.mkdir(parents=True, exist_ok=True)
        Path(support_output).parent.mkdir(parents=True, exist_ok=True)
        with bed_tmp.open("w", encoding="utf-8") as bed_handle, support_tmp.open(
            "w", encoding="utf-8"
        ) as support_handle:
            support_handle.write("\t".join(SUPPORT_COLUMNS) + "\n")
            current: tuple | None = None
            sample_counts: dict[str, int] = {}

            def emit(coords: tuple, observed_values: dict[str, int]) -> None:
                values = {
                    sample: int(observed_values.get(sample, 0))
                    for sample in expected_samples
                }
                supporting = sorted(
                    sample for sample, count in values.items() if count >= min_count
                )
                expected_support = candidate_support_by_key[coords]
                if len(supporting) != expected_support:
                    raise ValueError(
                        f"Candidate {coords} has {len(supporting)} supporting all-table "
                        f"samples but {expected_support} eligible-table samples"
                    )
                bed_handle.write("\t".join(map(str, coords)) + "\n")
                support_handle.write("\t".join([
                    *map(str, coords),
                    str(len(supporting)),
                    ",".join(supporting),
                    json.dumps(values, sort_keys=True, separators=(",", ":")),
                    "0",
                    "fithichip_all_interactions",
                    str(min_count),
                    str(min_samples),
                ]) + "\n")

            for row in conn.execute(query):
                coords = tuple(row[:6])
                if current is not None and coords != current:
                    emit(current, sample_counts)
                    sample_counts = {}
                current = coords
                if row[6] is not None:
                    sample_counts[str(row[6])] = int(row[7])
            if current is not None:
                emit(current, sample_counts)
        bed_tmp.replace(bedpe_output)
        support_tmp.replace(support_output)

        manifest = {
            "schema": "oracle-differential-hypothesis-universe-v1",
            "comparison": comparison,
            "source_kind": "fithichip_all_interactions",
            "fithichip_q_filter": None,
            "merge_nearby": False,
            "condition_labels_used_for_selection": False,
            "tolerance_bins": 0,
            "bin_size": bin_size,
            "lower_distance": lower_distance,
            "upper_distance": upper_distance,
            "min_count": min_count,
            "min_samples": min_samples,
            "samples": expected_samples,
            "eligible_rows_by_sample": eligible_rows,
            "all_interaction_rows_by_sample": all_rows,
            "candidate_count": candidate_count,
            "union_bedpe_sha256": _sha256(bedpe_output),
            "candidate_support_sha256": _sha256(support_output),
            "source_audits": {
                sample: {
                    "source_relative": audits[sample]["source_relative"],
                    "input_rows": audits[sample]["input_rows"],
                    "retained_all_rows": audits[sample]["retained_all_rows"],
                    "retained_abundance_eligible_rows": audits[sample][
                        "retained_abundance_eligible_rows"
                    ],
                }
                for sample in expected_samples
            },
        }
        temp_manifest = manifest_target.with_name(
            f".{manifest_target.name}.tmp-{os.getpid()}"
        )
        temp_manifest.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        temp_manifest.replace(manifest_target)
        log.info(
            "comparison=%s eligible_rows=%d all_rows_scanned=%d candidates=%d "
            "min_count=%d min_samples=%d",
            comparison, sum(eligible_rows.values()), sum(all_rows.values()),
            candidate_count, min_count, min_samples,
        )
        return manifest
    finally:
        if conn is not None:
            conn.close()
        database.unlink(missing_ok=True)
        bed_tmp.unlink(missing_ok=True)
        support_tmp.unlink(missing_ok=True)


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    build_universe(
        list(snakemake.input.eligible),
        list(snakemake.input.all_interactions),
        list(snakemake.input.audits),
        snakemake.input.blacklist,
        snakemake.input.chromsizes,
        snakemake.output.bedpe,
        snakemake.output.support,
        snakemake.output.manifest,
        comparison=str(snakemake.wildcards.comparison),
        expected_samples=list(snakemake.params.samples),
        bin_size=int(snakemake.params.bin_size),
        lower_distance=int(snakemake.params.lower_distance),
        upper_distance=int(snakemake.params.upper_distance),
        min_count=int(snakemake.params.min_count),
        min_samples=int(snakemake.params.min_samples),
    )


if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
