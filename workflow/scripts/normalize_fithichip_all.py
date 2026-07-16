"""Stream FitHiChIP's unmerged all-interaction table into a strict compact contract.

The official differential input is ``PREFIX.interactions_FitHiC.bed``.  It is
not the q-filtered or MergeNearContacts file.  Real libraries contain millions
of rows, so this module validates and writes chunks without loading the source
table into memory.  The ``eligible`` table applies only the predeclared raw
contact-count floor; no FitHiChIP p- or q-value influences eligibility.
"""
import gzip
import io
import json
import logging
import os
import sys
from contextlib import ExitStack, contextmanager
from pathlib import Path
from typing import Iterator, TextIO

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
try:
    from utils import open_text_auto, setup_logging  # type: ignore[attr-defined]
except ImportError:
    from utils import setup_logging  # type: ignore[no-redef]

    def open_text_auto(path: str | Path, mode: str = "rt") -> TextIO:
        """Temporary fallback until the shared signature-aware opener is available."""
        if mode != "rt":
            raise ValueError("fallback open_text_auto supports text reading only")
        source = Path(path)
        with source.open("rb") as handle:
            is_gzip = handle.read(2) == b"\x1f\x8b"
        if is_gzip:
            return gzip.open(source, mode, encoding="utf-8")
        return source.open(mode, encoding="utf-8")


log = logging.getLogger(__name__)

RAW_COORDS = ["chr1", "s1", "e1", "chr2", "s2", "e2"]
RAW_REQUIRED = RAW_COORDS + [
    "cc", "isPeak1", "isPeak2", "P-Value_Bias", "Q-Value_Bias",
]
OUT_COORDS = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
OUT_COLUMNS = OUT_COORDS + ["score", "pvalue", "fdr"]
PRIMARY = r"chr(?:[1-9]|1[0-9]|2[0-2]|X)"


@contextmanager
def _deterministic_gzip_text(path: str | Path) -> Iterator[TextIO]:
    """Atomically write deterministic gzip text (mtime=0)."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("wb") as raw:
            with gzip.GzipFile(
                filename="", fileobj=raw, mode="wb", mtime=0
            ) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                    yield text
        temporary.replace(target)
    finally:
        temporary.unlink(missing_ok=True)


def _read_chromsizes(path: str | Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with open_text_auto(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 2:
                raise ValueError(f"Malformed chromsizes row: {line.rstrip()!r}")
            sizes[fields[0]] = int(fields[1])
    if not sizes:
        raise ValueError("chromsizes is empty")
    return sizes


def _read_interval_index(
    path: str | Path,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    with open_text_auto(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                raise ValueError(f"Malformed blacklist row: {line.rstrip()!r}")
            chrom, start, end = fields[:3]
            start_i, end_i = int(start), int(end)
            if start_i < 0 or end_i <= start_i:
                raise ValueError(f"Invalid blacklist interval: {line.rstrip()!r}")
            by_chrom.setdefault(chrom, []).append((start_i, end_i))
    index: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for chrom, intervals in by_chrom.items():
        ordered = np.asarray(sorted(intervals), dtype=np.int64)
        index[chrom] = (ordered[:, 0], np.maximum.accumulate(ordered[:, 1]))
    return index


def _overlap_mask(
    chroms: pd.Series,
    starts: pd.Series,
    ends: pd.Series,
    index: dict[str, tuple[np.ndarray, np.ndarray]],
) -> np.ndarray:
    hit = np.zeros(len(chroms), dtype=bool)
    chrom_values = chroms.astype(str).to_numpy()
    start_values = starts.to_numpy(dtype=np.int64)
    end_values = ends.to_numpy(dtype=np.int64)
    for chrom in np.unique(chrom_values):
        if chrom not in index:
            continue
        rows = np.flatnonzero(chrom_values == chrom)
        bl_starts, prefix_max_end = index[chrom]
        pos = np.searchsorted(bl_starts, end_values[rows], side="left") - 1
        valid = pos >= 0
        hit[rows[valid]] = prefix_max_end[pos[valid]] > start_values[rows[valid]]
    return hit


def _validate_header(path: str | Path) -> list[str]:
    with open_text_auto(path, "rt") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                header = line.rstrip("\r\n").split("\t")
                break
        else:
            raise ValueError("FitHiChIP all-interaction table is empty")
    if header[:6] != RAW_COORDS:
        raise ValueError(
            "Expected FitHiChIP all-interaction header chr1/s1/e1/chr2/s2/e2; "
            f"found {header[:6]}"
        )
    missing = [column for column in RAW_REQUIRED if column not in header]
    if missing:
        raise ValueError(f"FitHiChIP all-interaction table lacks columns: {missing}")
    if len(header) != len(set(header)):
        raise ValueError("FitHiChIP all-interaction header contains duplicate columns")
    return header


def _strict_numeric(chunk: pd.DataFrame, columns: list[str]) -> None:
    for column in columns:
        converted = pd.to_numeric(chunk[column], errors="coerce")
        if converted.isna().any():
            bad = int(converted.isna().sum())
            raise ValueError(f"FitHiChIP column {column!r} has {bad} non-numeric values")
        chunk[column] = converted


def _strict_probability_or_missing(values: pd.Series, column: str) -> pd.Series:
    """Accept numeric probabilities or FitHiChIP/R's explicit NA tokens only."""
    raw = values.astype(str).str.strip()
    missing = raw.str.casefold().isin({"na", "nan"})
    numeric = pd.to_numeric(raw.where(~missing), errors="coerce")
    malformed = (~missing) & numeric.isna()
    out_of_range = (~missing) & (
        (~np.isfinite(numeric)) | (numeric < 0) | (numeric > 1)
    )
    if malformed.any() or out_of_range.any():
        examples = raw.loc[malformed | out_of_range].drop_duplicates().head(5).tolist()
        raise ValueError(
            f"FitHiChIP column {column!r} must contain probabilities in [0,1] "
            f"or explicit NA/NaN tokens; invalid values: {examples}"
        )
    return numeric.astype(float)


def _validate_native_geometry(
    chunk: pd.DataFrame,
    chromsizes: dict[str, int],
    bin_size: int,
    lower_distance: int,
    upper_distance: int,
    interaction_type: str,
) -> None:
    numeric = ["s1", "e1", "s2", "e2", "cc", "isPeak1", "isPeak2"]
    _strict_numeric(chunk, numeric)
    integer = chunk[numeric].apply(lambda col: np.equal(col, np.floor(col)))
    if not bool(integer.to_numpy().all()):
        raise ValueError("FitHiChIP coordinates, cc, and peak flags must be integers")
    chunk[numeric] = chunk[numeric].astype(np.int64)

    if (chunk["cc"] <= 0).any():
        raise ValueError("FitHiChIP all-interaction cc values must be positive")
    if not chunk["isPeak1"].isin([0, 1]).all() or not chunk["isPeak2"].isin([0, 1]).all():
        raise ValueError("FitHiChIP peak flags must be 0 or 1")

    bad_geometry = (
        (chunk["s1"] < 0)
        | (chunk["s2"] < 0)
        | (chunk["e1"] - chunk["s1"] != bin_size)
        | (chunk["e2"] - chunk["s2"] != bin_size)
        | (chunk["s1"] % bin_size != 0)
        | (chunk["s2"] % bin_size != 0)
        | (chunk["chr1"].astype(str) != chunk["chr2"].astype(str))
        | (chunk["s1"] >= chunk["s2"])
    )
    distance = chunk["s2"] - chunk["s1"]
    bad_geometry |= (distance < lower_distance) | (distance > upper_distance)
    chrom_names = chunk["chr1"].astype(str)
    unknown = ~chrom_names.isin(chromsizes)
    unknown_primary = unknown & chrom_names.str.fullmatch(PRIMARY)
    if unknown_primary.any():
        names = sorted(chunk.loc[unknown_primary, "chr1"].astype(str).unique())
        raise ValueError(
            f"FitHiChIP primary chromosomes are absent from chromsizes: {names}"
        )
    known = ~unknown
    if known.any():
        chrom_end = chunk.loc[known, "chr1"].map(chromsizes).astype(np.int64)
        bad_geometry.loc[known] |= (
            (chunk.loc[known, "e1"] > chrom_end)
            | (chunk.loc[known, "e2"] > chrom_end)
        )
    if bad_geometry.any():
        raise ValueError(
            f"FitHiChIP all-interaction table has {int(bad_geometry.sum())} rows "
            "outside its native bin/cis/distance/chromosome contract"
        )

    peak1, peak2 = chunk["isPeak1"], chunk["isPeak2"]
    if interaction_type == "Peak-to-Peak":
        expected = (peak1 == 1) & (peak2 == 1)
    elif interaction_type == "Peak-to-NonPeak":
        expected = peak1 != peak2
    elif interaction_type == "Peak-to-ALL":
        expected = (peak1 == 1) | (peak2 == 1)
    elif interaction_type == "ALL-to-ALL":
        expected = pd.Series(True, index=chunk.index)
    else:
        raise ValueError(f"Unsupported FitHiChIP interaction type: {interaction_type!r}")
    if not expected.all():
        raise ValueError(
            f"FitHiChIP all-interaction rows violate {interaction_type} peak flags"
        )


def normalize_all_interactions(
    source: str | Path,
    blacklist: str | Path,
    chromsizes_path: str | Path,
    all_output: str | Path,
    eligible_output: str | Path,
    audit_output: str | Path,
    *,
    sample: str,
    bin_size: int,
    lower_distance: int,
    upper_distance: int,
    min_count: int,
    interaction_type: str,
    source_relative: str,
    chunk_size: int = 250_000,
) -> dict:
    """Normalize a native FitHiChIP table and return its machine-readable audit."""
    if min_count < 1:
        raise ValueError("differential min_count must be at least one")
    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")
    header = _validate_header(source)
    chromsizes = _read_chromsizes(chromsizes_path)
    blacklist_index = _read_interval_index(blacklist)

    counts = {
        "input_rows": 0,
        "removed_non_primary": 0,
        "removed_blacklist": 0,
        "retained_all_rows": 0,
        "retained_abundance_eligible_rows": 0,
        "retained_q_gt_0_05_rows": 0,
        "retained_missing_q_rows": 0,
    }
    first_all = True
    first_eligible = True
    with ExitStack() as stack:
        all_handle = stack.enter_context(_deterministic_gzip_text(all_output))
        eligible_handle = stack.enter_context(_deterministic_gzip_text(eligible_output))
        reader = pd.read_csv(
            source,
            sep="\t",
            usecols=RAW_REQUIRED,
            chunksize=chunk_size,
            dtype={
                "chr1": str,
                "chr2": str,
                "P-Value_Bias": str,
                "Q-Value_Bias": str,
            },
            keep_default_na=False,
        )
        for chunk in reader:
            counts["input_rows"] += len(chunk)
            _validate_native_geometry(
                chunk, chromsizes, bin_size, lower_distance, upper_distance,
                interaction_type,
            )
            for column in ("P-Value_Bias", "Q-Value_Bias"):
                chunk[column] = _strict_probability_or_missing(chunk[column], column)

            primary = chunk["chr1"].astype(str).str.fullmatch(PRIMARY)
            counts["removed_non_primary"] += int((~primary).sum())
            chunk = chunk.loc[primary].copy()

            blocked = _overlap_mask(
                chunk["chr1"], chunk["s1"], chunk["e1"], blacklist_index
            ) | _overlap_mask(
                chunk["chr2"], chunk["s2"], chunk["e2"], blacklist_index
            )
            counts["removed_blacklist"] += int(blocked.sum())
            chunk = chunk.loc[~blocked].copy()

            out = chunk.rename(columns={
                "chr1": "chrom1", "s1": "start1", "e1": "end1",
                "chr2": "chrom2", "s2": "start2", "e2": "end2",
                "cc": "score", "P-Value_Bias": "pvalue",
                "Q-Value_Bias": "fdr",
            })[OUT_COLUMNS]
            counts["retained_all_rows"] += len(out)
            counts["retained_missing_q_rows"] += int(out["fdr"].isna().sum())
            counts["retained_q_gt_0_05_rows"] += int(
                (out["fdr"].notna() & (out["fdr"] > 0.05)).sum()
            )
            out.to_csv(all_handle, sep="\t", index=False, header=first_all)
            first_all = False

            eligible = out.loc[out["score"] >= min_count]
            counts["retained_abundance_eligible_rows"] += len(eligible)
            eligible.to_csv(
                eligible_handle, sep="\t", index=False, header=first_eligible
            )
            first_eligible = False

        if first_all:
            pd.DataFrame(columns=OUT_COLUMNS).to_csv(
                all_handle, sep="\t", index=False
            )
        if first_eligible:
            pd.DataFrame(columns=OUT_COLUMNS).to_csv(
                eligible_handle, sep="\t", index=False
            )

    if counts["input_rows"] == 0 or counts["retained_all_rows"] == 0:
        raise RuntimeError("No primary, blacklist-clean FitHiChIP interactions remained")
    if counts["input_rows"] != (
        counts["removed_non_primary"]
        + counts["removed_blacklist"]
        + counts["retained_all_rows"]
    ):
        raise AssertionError("FitHiChIP normalization audit does not conserve rows")

    audit = {
        "schema": "oracle-fithichip-all-interactions-v1",
        "sample": sample,
        "source_kind": "fithichip_all_interactions",
        "source_relative": source_relative,
        "source_header": header,
        "fithichip_q_filter": None,
        "merge_nearby": False,
        "contact_count_filter_all_table": None,
        "eligible_min_count": min_count,
        "bin_size": bin_size,
        "lower_distance": lower_distance,
        "upper_distance": upper_distance,
        "interaction_type": interaction_type,
        **counts,
    }
    target = Path(audit_output)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.tmp-{os.getpid()}")
    try:
        temporary.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
        temporary.replace(target)
    finally:
        temporary.unlink(missing_ok=True)
    log.info(
        "sample=%s input=%d retained_all=%d eligible_cc_ge_%d=%d",
        sample, counts["input_rows"], counts["retained_all_rows"], min_count,
        counts["retained_abundance_eligible_rows"],
    )
    return audit


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    normalize_all_interactions(
        snakemake.input.raw,
        snakemake.input.blacklist,
        snakemake.input.chromsizes,
        snakemake.output.all_interactions,
        snakemake.output.eligible,
        snakemake.output.audit,
        sample=str(snakemake.wildcards.sample),
        bin_size=int(snakemake.params.bin_size),
        lower_distance=int(snakemake.params.lower_distance),
        upper_distance=int(snakemake.params.upper_distance),
        min_count=int(snakemake.params.min_count),
        interaction_type=str(snakemake.params.interaction_type),
        source_relative=str(snakemake.params.source_relative),
    )


if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
