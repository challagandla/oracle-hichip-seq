"""Paired differential testing on a predeclared FitHiChIP pixel universe.

The hypothesis universe is built upstream from unmerged, q-unthresholded
``PREFIX.interactions_FitHiC.bed`` tables.  This script revalidates that contract,
requires every cooler-derived count table to contain the exact same ordered
pixels, and then fits ``~ pairing_factor + condition`` with pyDESeq2.  Small
paired cohorts may run as explicitly labelled pilot analyses, but cannot be
mistaken for publication-ready inference.
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils import setup_logging  # noqa: E402


LOOP_COORD_COLS = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
SUPPORT_COLUMNS = LOOP_COORD_COLS + [
    "sample_support", "support_samples", "source_counts", "tolerance_bins",
    "source_kind", "min_count", "min_samples",
]


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _loop_keys(coords: pd.DataFrame) -> pd.Series:
    encoded = coords[LOOP_COORD_COLS].astype(str).agg("\x1f".join, axis=1)
    return encoded.map(lambda value: hashlib.sha1(value.encode("utf-8")).hexdigest())


def _strict_coordinates(frame: pd.DataFrame, label: str) -> pd.DataFrame:
    missing = [column for column in LOOP_COORD_COLS if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} lacks coordinate columns: {missing}")
    coords = frame[LOOP_COORD_COLS].copy()
    if coords[["chrom1", "chrom2"]].isna().any().any():
        raise ValueError(f"{label} has missing chromosome values")
    coords["chrom1"] = coords["chrom1"].astype(str)
    coords["chrom2"] = coords["chrom2"].astype(str)
    for column in ("start1", "end1", "start2", "end2"):
        numeric = pd.to_numeric(coords[column], errors="coerce")
        if numeric.isna().any() or not np.equal(numeric, np.floor(numeric)).all():
            raise ValueError(f"{label} has non-integer {column} values")
        coords[column] = numeric.astype(np.int64)
    bad = (
        (coords["chrom1"] != coords["chrom2"])
        | (coords["start1"] < 0)
        | (coords["start2"] < 0)
        | (coords["end1"] <= coords["start1"])
        | (coords["end2"] <= coords["start2"])
        | (coords["start1"] >= coords["start2"])
    )
    if bad.any():
        raise ValueError(f"{label} has {int(bad.sum())} invalid cis upper-triangle rows")
    return coords


def _load_universe(path: str | Path) -> pd.DataFrame:
    try:
        raw = pd.read_csv(path, sep="\t", header=None, dtype=str)
    except pd.errors.EmptyDataError as exc:
        raise ValueError("Differential hypothesis universe is empty") from exc
    if raw.shape[1] != 6:
        raise ValueError(
            f"Differential hypothesis universe must be headerless six-column BEDPE; "
            f"found {raw.shape[1]} columns"
        )
    raw.columns = LOOP_COORD_COLS
    coords = _strict_coordinates(raw, "differential hypothesis universe")
    keys = _loop_keys(coords)
    if keys.duplicated().any():
        raise ValueError("Differential hypothesis universe contains duplicate pixels")
    coords.insert(0, "loop_key", keys)
    return coords


def _coordinates_equal(left: pd.DataFrame, right: pd.DataFrame) -> bool:
    return (
        len(left) == len(right)
        and left[LOOP_COORD_COLS].reset_index(drop=True).equals(
            right[LOOP_COORD_COLS].reset_index(drop=True)
        )
    )


def _load_count_table(
    files: list[str],
    expected_coords: pd.DataFrame | None = None,
    expected_samples: list[str] | None = None,
) -> tuple[pd.DataFrame, list[str], pd.DataFrame]:
    """Load exact, complete integer count vectors for every selected sample."""
    if not files:
        raise ValueError("No differential count tables were supplied")
    if expected_samples is not None and len(files) != len(expected_samples):
        raise ValueError("Count-table count must equal the selected sample count")

    vectors: list[pd.Series] = []
    sample_order: list[str] = []
    reference = expected_coords.copy() if expected_coords is not None else None
    for position, filename in enumerate(files):
        frame = pd.read_csv(filename, sep="\t")
        required = LOOP_COORD_COLS + ["count", "sample"]
        missing = [column for column in required if column not in frame.columns]
        if missing:
            raise ValueError(f"Count table {filename} lacks columns: {missing}")
        if frame.empty:
            raise ValueError(f"Count table {filename} contains no loop rows")
        if frame["sample"].isna().any() or (
            frame["sample"].astype(str).str.strip() == ""
        ).any():
            raise ValueError(f"Count table {filename} has missing sample IDs")
        sample_values = frame["sample"].astype(str).unique().tolist()
        if len(sample_values) != 1 or not sample_values[0]:
            raise ValueError(f"Count table {filename} must contain exactly one sample ID")
        sample = sample_values[0]
        if not (frame["sample"].astype(str) == sample).all():
            raise ValueError(f"Count table {filename} contains mixed sample IDs")
        if sample in sample_order:
            raise ValueError(f"Duplicate differential count table for sample {sample!r}")
        if expected_samples is not None:
            expected_sample = str(expected_samples[position])
            if sample != expected_sample:
                raise ValueError(
                    f"Count table {filename} contains sample {sample!r}; expected "
                    f"{expected_sample!r} at this position"
                )
            if not Path(filename).name.startswith(f"{expected_sample}."):
                raise ValueError(
                    f"Count-table filename {filename} does not match sample "
                    f"{expected_sample!r}"
                )

        coords = _strict_coordinates(frame, f"count table for {sample}")
        keys = _loop_keys(coords)
        if keys.duplicated().any():
            raise ValueError(f"Count table for {sample} contains duplicate pixels")
        if reference is None:
            reference = coords.copy()
            reference.insert(0, "loop_key", keys)
        elif not _coordinates_equal(coords, reference):
            raise ValueError(
                f"Count table for {sample} does not contain the exact ordered "
                "hypothesis universe"
            )

        values = pd.to_numeric(frame["count"], errors="coerce")
        if (
            values.isna().any()
            or not np.isfinite(values).all()
            or not np.equal(values, np.floor(values)).all()
            or (values < 0).any()
        ):
            raise ValueError(
                f"Count table for {sample} must contain finite non-negative integers"
            )
        vectors.append(
            pd.Series(values.astype(np.int64).to_numpy(), index=keys, name=sample)
        )
        sample_order.append(sample)

    if reference is None:
        raise RuntimeError("No count-table coordinates were loaded")
    matrix = pd.concat(vectors, axis=1)
    if matrix.isna().any().any() or len(matrix) != len(reference):
        raise AssertionError("Count matrix lost its exact hypothesis-universe alignment")
    return matrix, sample_order, reference


def _load_source_audits(
    paths: list[str], expected_samples: list[str], min_count: int
) -> dict[str, dict]:
    audits: dict[str, dict] = {}
    for path in paths:
        audit = json.loads(Path(path).read_text(encoding="utf-8"))
        sample = str(audit.get("sample", ""))
        if not sample or sample in audits:
            raise ValueError(f"Duplicate or empty sample in source audit {path}")
        expected = {
            "schema": "oracle-fithichip-all-interactions-v1",
            "source_kind": "fithichip_all_interactions",
            "fithichip_q_filter": None,
            "merge_nearby": False,
            "eligible_min_count": min_count,
        }
        for key, value in expected.items():
            if audit.get(key) != value:
                raise ValueError(
                    f"Source audit for {sample} has {key}={audit.get(key)!r}; "
                    f"expected {value!r}"
                )
        audits[sample] = audit
    if set(audits) != set(expected_samples):
        raise ValueError(
            "Source-audit samples do not equal selected samples; "
            f"expected={sorted(expected_samples)}, observed={sorted(audits)}"
        )
    return audits


def _parse_source_counts(value: object, row_number: int) -> dict[str, int]:
    try:
        raw = json.loads(str(value))
    except json.JSONDecodeError as exc:
        raise ValueError(f"candidate_support row {row_number} has invalid source_counts") from exc
    if not isinstance(raw, dict):
        raise ValueError(f"candidate_support row {row_number} source_counts is not an object")
    parsed: dict[str, int] = {}
    for sample, count in raw.items():
        if (
            isinstance(count, bool)
            or not isinstance(count, int)
            or count < 0
            or not str(sample)
        ):
            raise ValueError(
                f"candidate_support row {row_number} has invalid source count "
                f"{sample!r}: {count!r}"
            )
        parsed[str(sample)] = count
    return parsed


def validate_hypothesis_contract(
    universe_path: str | Path,
    support_path: str | Path,
    manifest_path: str | Path,
    audit_paths: list[str],
    expected_samples: list[str],
    min_count: int,
    min_samples: int,
    matrix: pd.DataFrame,
    coords: pd.DataFrame,
) -> dict:
    """Prove that candidate selection was exact-grid, label-blind, and q-free."""
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    expected_manifest = {
        "schema": "oracle-differential-hypothesis-universe-v1",
        "source_kind": "fithichip_all_interactions",
        "fithichip_q_filter": None,
        "merge_nearby": False,
        "condition_labels_used_for_selection": False,
        "tolerance_bins": 0,
        "min_count": min_count,
        "min_samples": min_samples,
        "candidate_count": len(coords),
        "samples": expected_samples,
        "union_bedpe_sha256": _sha256(universe_path),
        "candidate_support_sha256": _sha256(support_path),
    }
    for key, value in expected_manifest.items():
        if manifest.get(key) != value:
            raise ValueError(
                f"Hypothesis-universe manifest has {key}={manifest.get(key)!r}; "
                f"expected {value!r}"
            )
    _load_source_audits(audit_paths, expected_samples, min_count)

    support = pd.read_csv(support_path, sep="\t", dtype=str, keep_default_na=False)
    missing = [column for column in SUPPORT_COLUMNS if column not in support.columns]
    if missing:
        raise ValueError(f"candidate_support table lacks columns: {missing}")
    support_coords = _strict_coordinates(support, "candidate_support table")
    if not _coordinates_equal(support_coords, coords):
        raise ValueError("candidate_support rows do not equal the ordered hypothesis universe")

    for row_number, (row, key) in enumerate(
        zip(support.itertuples(index=False), coords["loop_key"]), start=2
    ):
        record = row._asdict()
        try:
            sample_support = int(record["sample_support"])
            tolerance = int(record["tolerance_bins"])
            row_min_count = int(record["min_count"])
            row_min_samples = int(record["min_samples"])
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"candidate_support row {row_number} has non-integer contract fields"
            ) from exc
        support_samples = [
            sample for sample in str(record["support_samples"]).split(",") if sample
        ]
        source_counts = _parse_source_counts(record["source_counts"], row_number)
        derived_support = {
            sample for sample, count in source_counts.items() if count >= min_count
        }
        if (
            record["source_kind"] != "fithichip_all_interactions"
            or tolerance != 0
            or row_min_count != min_count
            or row_min_samples != min_samples
            or len(support_samples) != len(set(support_samples))
            or set(source_counts) != set(expected_samples)
            or set(support_samples) != derived_support
            or sample_support != len(derived_support)
            or sample_support < min_samples
        ):
            raise ValueError(f"candidate_support row {row_number} violates its q-free contract")
        observed = {
            sample: int(matrix.loc[key, sample])
            for sample in expected_samples
        }
        if source_counts != observed:
            raise ValueError(
                f"candidate_support row {row_number} disagrees with cooler counts; "
                f"source={source_counts}, cooler={observed}"
            )
    return manifest


def build_design_metadata(
    columns: list[str], cases: list[str], controls: list[str], samples: pd.DataFrame,
    paired_by: str | None, covariates: list[str],
) -> tuple[pd.DataFrame, list[str]]:
    """Construct and validate a full-rank paired differential design."""
    if set(cases) & set(controls):
        raise ValueError("Case and control sample sets overlap")
    if len(columns) != len(set(columns)) or set(columns) != set(cases + controls):
        raise ValueError("Count matrix columns must equal the configured case/control samples")
    if not samples.index.is_unique:
        raise ValueError("samples.tsv contains duplicate sample_id values")
    missing = sorted(set(columns) - set(samples.index))
    if missing:
        raise ValueError(f"samples.tsv lacks differential samples: {missing}")

    metadata = samples.loc[columns].copy()
    metadata["condition"] = ["case" if sid in cases else "control" for sid in columns]
    design_factors: list[str] = []
    if paired_by:
        if paired_by not in metadata.columns:
            raise ValueError(f"paired_by column {paired_by!r} is absent from samples.tsv")
        if metadata[paired_by].isna().any() or (
            metadata[paired_by].astype(str).str.strip() == ""
        ).any():
            raise ValueError(f"paired_by column {paired_by!r} contains empty values")
        pairing = metadata.groupby(paired_by)["condition"].value_counts().unstack(fill_value=0)
        bad = pairing[(pairing.get("case", 0) != 1) | (pairing.get("control", 0) != 1)]
        if not bad.empty:
            raise ValueError(
                f"Comparison is not one-to-one paired by {paired_by}; unmatched levels: "
                f"{bad.index.astype(str).tolist()}"
            )
        if len(pairing) < 2:
            raise ValueError("Paired differential analysis requires at least two complete pairs")
        design_factors.append(paired_by)

    for factor in covariates:
        if factor not in metadata.columns:
            raise ValueError(f"Differential covariate {factor!r} is absent from samples.tsv")
        if metadata[factor].nunique(dropna=False) > 1:
            design_factors.append(factor)
    design_factors.append("condition")
    design_factors = list(dict.fromkeys(design_factors))
    design_metadata = metadata[design_factors]
    matrix = pd.get_dummies(design_metadata.astype(str), drop_first=True, dtype=float)
    matrix.insert(0, "intercept", 1.0)
    if np.linalg.matrix_rank(matrix.to_numpy()) < matrix.shape[1]:
        raise ValueError(
            "Differential design is rank-deficient; a covariate is confounded "
            "with condition or the pairing factor"
        )
    return design_metadata, design_factors


def _pairing_summary(
    metadata: pd.DataFrame, paired_by: str, cases: list[str], controls: list[str]
) -> list[dict[str, str]]:
    pairs: list[dict[str, str]] = []
    for pair_id, group in metadata.groupby(paired_by, sort=True):
        case_ids = [sample for sample in group.index if sample in cases]
        control_ids = [sample for sample in group.index if sample in controls]
        if len(case_ids) != 1 or len(control_ids) != 1:
            raise ValueError(f"Pair {pair_id!r} is not one-to-one case/control")
        pairs.append({
            "pair_id": str(pair_id),
            "case_sample": case_ids[0],
            "control_sample": control_ids[0],
        })
    return pairs


def classify_analysis_status(
    n_complete_pairs: int,
    publication_min_complete_pairs: int,
    require_publication_ready: bool,
) -> tuple[str, bool]:
    """Classify a contrast and optionally enforce the hard publication gate."""
    if n_complete_pairs < 2:
        raise ValueError("Paired differential analysis requires at least two complete pairs")
    if publication_min_complete_pairs < 3:
        raise ValueError("publication minimum must be at least three complete pairs")
    publication_eligible = n_complete_pairs >= publication_min_complete_pairs
    if require_publication_ready and not publication_eligible:
        raise RuntimeError(
            f"Comparison has {n_complete_pairs} complete pairs, below the configured "
            f"publication minimum {publication_min_complete_pairs}"
        )
    return (
        "STANDARD_INFERENCE" if publication_eligible else "PILOT_UNDERPOWERED",
        publication_eligible,
    )


def _pydeseq2(
    matrix: pd.DataFrame, metadata: pd.DataFrame, design_factors: list[str],
    fdr: float, lfc: float, n_cpus: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    if metadata["condition"].nunique() != 2:
        raise RuntimeError("Differential analysis requires both case and control samples")
    counts = matrix.T
    dds = DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design_factors=design_factors,
        ref_level=["condition", "control"],
        refit_cooks=False,
        n_cpus=n_cpus,
        quiet=True,
    )
    dds.deseq2()
    stat = DeseqStats(dds, contrast=("condition", "case", "control"), quiet=True)
    stat.summary()
    result = stat.results_df.copy()
    result["sig"] = (
        (result["padj"] < fdr) & (result["log2FoldChange"].abs() >= lfc)
    ).fillna(False)
    normalized_array = dds.layers.get("normed_counts")
    if normalized_array is None:
        raise RuntimeError("pyDESeq2 did not expose its normalized-count layer")
    normalized = pd.DataFrame(
        np.asarray(normalized_array, dtype=float),
        index=counts.index.astype(str),
        columns=counts.columns.astype(str),
    ).T
    return result, normalized


def build_paired_effects(
    normalized: pd.DataFrame,
    coords: pd.DataFrame,
    pairs: list[dict[str, str]],
    paired_by: str,
    analysis_status: str,
) -> pd.DataFrame:
    """Return per-loop, per-pair normalized effects for audit and forest plots."""
    if set(normalized.index.astype(str)) != set(coords["loop_key"].astype(str)):
        raise ValueError("Normalized-count rows do not match the hypothesis universe")
    normalized = normalized.loc[coords["loop_key"].astype(str)]
    frames: list[pd.DataFrame] = []
    for pair in pairs:
        case_sample = pair["case_sample"]
        control_sample = pair["control_sample"]
        if case_sample not in normalized.columns or control_sample not in normalized.columns:
            raise ValueError(f"Normalized counts are missing pair {pair['pair_id']!r}")
        case_values = normalized[case_sample].to_numpy(dtype=float)
        control_values = normalized[control_sample].to_numpy(dtype=float)
        frame = coords.copy()
        frame["pairing_factor"] = paired_by
        frame["pair_id"] = pair["pair_id"]
        frame[paired_by] = pair["pair_id"]
        frame["case_sample"] = case_sample
        frame["control_sample"] = control_sample
        frame["case_normalized_count"] = case_values
        frame["control_normalized_count"] = control_values
        frame["paired_log2_ratio"] = np.log2(
            (case_values + 1.0) / (control_values + 1.0)
        )
        frame["analysis_status"] = analysis_status
        frames.append(frame)
    if not frames:
        raise ValueError("No complete pairs were available for paired-effect output")
    return pd.concat(frames, ignore_index=True)


def _volcano(result: pd.DataFrame, out_png: str | Path, fdr: float, lfc: float) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5, 4.5))
    x = result["log2FoldChange"]
    y = -np.log10(result["padj"].replace(0, np.nextafter(0, 1)))
    sig = result["sig"].astype(bool)
    ax.scatter(x[~sig], y[~sig], s=4, c="#ccc", alpha=0.5)
    ax.scatter(x[sig & (x > 0)], y[sig & (x > 0)], s=6, c="#c0392b", label="increased")
    ax.scatter(x[sig & (x < 0)], y[sig & (x < 0)], s=6, c="#1f5fbf", label="decreased")
    ax.axhline(-np.log10(fdr), c="k", ls="--", lw=0.5)
    ax.axvline(lfc, c="k", ls="--", lw=0.5)
    ax.axvline(-lfc, c="k", ls="--", lw=0.5)
    ax.set_xlabel("log2 fold change (case / control)")
    ax.set_ylabel("-log10 adjusted p")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def _ma(result: pd.DataFrame, out_png: str | Path) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5, 4))
    base = np.log10(result["baseMean"].replace(0, np.nextafter(0, 1)))
    ax.scatter(
        base, result["log2FoldChange"], s=4,
        c=np.where(result["sig"].astype(bool), "#c0392b", "#bbb"), alpha=0.6,
    )
    ax.axhline(0, c="k", lw=0.5)
    ax.set_xlabel("log10(base mean)")
    ax.set_ylabel("log2 fold change")
    fig.tight_layout()
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    cases_raw, controls_raw, comparison_config = snakemake.params.groups
    cases, controls = list(cases_raw), list(controls_raw)
    comparison_config = dict(comparison_config)
    expected_samples = cases + controls
    if not cases or not controls:
        raise RuntimeError(
            f"Comparison {snakemake.wildcards.comparison} has empty case/control groups"
        )

    paired_by = str(snakemake.params.paired_by) if snakemake.params.paired_by else ""
    if not paired_by:
        raise ValueError("Differential analysis requires an explicit paired_by column")
    covariates = list(snakemake.params.covariates)
    min_count = int(snakemake.params.min_count)
    min_samples = int(snakemake.params.min_samples)
    publication_min = int(snakemake.params.publication_min_complete_pairs)
    require_publication = bool(snakemake.params.require_publication_ready)

    coords = _load_universe(snakemake.input.universe)
    matrix, order, count_coords = _load_count_table(
        list(snakemake.input.counts), coords, expected_samples
    )
    if order != expected_samples or not _coordinates_equal(count_coords, coords):
        raise AssertionError("Differential counts lost configured sample/pixel order")
    universe_manifest = validate_hypothesis_contract(
        snakemake.input.universe,
        snakemake.input.support,
        snakemake.input.universe_manifest,
        list(snakemake.input.source_audits),
        expected_samples,
        min_count,
        min_samples,
        matrix,
        coords,
    )

    sample_sheet = pd.read_csv(
        snakemake.input.samples,
        sep="\t",
        comment="#",
        dtype=str,
        keep_default_na=False,
    ).set_index("sample_id")
    design_metadata, design_factors = build_design_metadata(
        list(matrix.columns), cases, controls, sample_sheet, paired_by, covariates
    )
    selected_metadata = sample_sheet.loc[list(matrix.columns)].copy()
    selected_metadata["condition"] = [
        "case" if sample in cases else "control" for sample in selected_metadata.index
    ]
    pairs = _pairing_summary(selected_metadata, paired_by, cases, controls)
    n_complete_pairs = len(pairs)
    analysis_status, publication_eligible = classify_analysis_status(
        n_complete_pairs, publication_min, require_publication
    )

    method = str(snakemake.params.method)
    if method != "pyDESeq2":
        raise ValueError(f"Unsupported differential method {method!r}; use pyDESeq2")
    fdr = float(snakemake.params.fdr)
    lfc = float(snakemake.params.log2fc_min)
    result, normalized = _pydeseq2(
        matrix,
        design_metadata,
        design_factors,
        fdr=fdr,
        lfc=lfc,
        n_cpus=int(getattr(snakemake, "threads", 1)),
    )
    result.index = result.index.astype(str)
    if set(result.index) != set(coords["loop_key"]):
        raise RuntimeError("pyDESeq2 results do not match the tested hypothesis universe")
    result.index.name = "loop_key"
    result = (
        result.reset_index()
        .merge(coords, on="loop_key", how="left", validate="one_to_one")
        .sort_values("padj", na_position="last")
    )
    result["analysis_status"] = analysis_status
    required_result = {"baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj", "sig"}
    missing_result = sorted(required_result - set(result.columns))
    if missing_result:
        raise RuntimeError(f"pyDESeq2 result lacks required columns: {missing_result}")

    paired_effects = build_paired_effects(
        normalized, coords, pairs, paired_by, analysis_status
    )
    Path(snakemake.output.tsv).parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(snakemake.output.tsv, sep="\t", index=False)
    paired_effects.to_csv(snakemake.output.paired_effects, sep="\t", index=False)

    design = {
        "comparison": str(snakemake.wildcards.comparison),
        "analysis_status": analysis_status,
        "n_complete_pairs": n_complete_pairs,
        "publication_eligible": publication_eligible,
        "publication_min_complete_pairs": publication_min,
        "require_publication_ready": require_publication,
        "paired_subjects": [pair["pair_id"] for pair in pairs],
        "pairing_factor": paired_by,
        "pairs": pairs,
        "cases": cases,
        "controls": controls,
        "covariates_requested": covariates,
        "design_factors_fitted": design_factors,
        "design_formula": "~ " + " + ".join(design_factors),
        "metadata": design_metadata.reset_index().rename(
            columns={"index": "sample_id"}
        ).to_dict("records"),
        "candidate_loops": int(len(matrix)),
        "tested_loops": int(len(matrix)),
        "prefilter": {
            "source": "fithichip_all_interactions",
            "q_value_filter": None,
            "merge_nearby": False,
            "tolerance_bins": 0,
            "min_count": min_count,
            "min_samples": min_samples,
        },
        "hypothesis_universe": universe_manifest,
        "config": comparison_config,
        "interpretation": (
            f"Differential {comparison_config.get('mark', 'mark')}-associated contact "
            "signal; changes can reflect anchor occupancy, contact frequency, or both."
        ),
    }
    Path(snakemake.output.design).write_text(
        json.dumps(design, indent=2) + "\n", encoding="utf-8"
    )
    _volcano(result, snakemake.output.volcano, fdr=fdr, lfc=lfc)
    _ma(result, snakemake.output.ma)


if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
