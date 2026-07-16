"""Fail-fast validation for metadata and cross-stage resolution/design contracts."""
import re
from pathlib import Path

import numpy as np
import pandas as pd

REQUIRED_SAMPLE_COLUMNS = {
    "sample_id", "subject_id", "cell_type", "tissue", "replicate", "mark",
    "fastq_r1", "fastq_r2", "batch", "library_protocol", "restriction_enzyme",
}
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
SAFE_PATH = re.compile(r"^[A-Za-z0-9_./:+-]+$")
SRA_ACCESSION = re.compile(r"^(?:SRR|ERR|DRR)[0-9]+$")


def _selected(
    samples: pd.DataFrame,
    filters: dict,
    mark: str | None,
    include_subjects: list[str] | None = None,
) -> pd.DataFrame:
    selected = samples
    if mark is not None:
        selected = selected[selected["mark"] == str(mark)]
    for column, value in (filters or {}).items():
        if column not in samples.columns:
            raise ValueError(f"comparison filter column {column!r} is absent from samples.tsv")
        selected = selected[selected[column] == str(value)]
    if include_subjects is not None:
        selected = selected[selected["subject_id"].isin(include_subjects)]
    return selected


def validate_pipeline_config(config: dict, genome_cfg: dict, samples: pd.DataFrame) -> None:
    errors: list[str] = []
    missing = sorted(REQUIRED_SAMPLE_COLUMNS - set(samples.columns))
    if missing:
        errors.append(f"samples.tsv missing required columns: {missing}")
    if errors:
        raise ValueError("Configuration validation failed:\n- " + "\n- ".join(errors))

    text = samples.fillna("").astype(str)
    group_by_raw = config.get("hicrep", {}).get("group_by", [])
    if not isinstance(group_by_raw, list) or not group_by_raw:
        errors.append("hicrep.group_by must be a non-empty list of sample-sheet columns")
        hicrep_group_by: list[str] = [
            column for column in ("cell_type", "mark") if column in text.columns
        ]
    else:
        hicrep_group_by = [str(column) for column in group_by_raw]
        if len(hicrep_group_by) != len(set(hicrep_group_by)):
            errors.append("hicrep.group_by must not contain duplicate columns")
        missing_group_columns = sorted(set(hicrep_group_by) - set(text.columns))
        if missing_group_columns:
            errors.append(
                f"hicrep.group_by columns are absent from samples.tsv: {missing_group_columns}"
            )
        forbidden_group_columns = sorted(
            set(hicrep_group_by)
            & {
                "sample_id", "subject_id", "replicate", "batch", "srr",
                "fastq_r1", "fastq_r2", "notes",
            }
        )
        if forbidden_group_columns:
            errors.append(
                "hicrep.group_by must describe biological condition/assay strata, "
                f"not identifiers or technical fields: {forbidden_group_columns}"
            )
        if "mark" not in hicrep_group_by:
            errors.append("hicrep.group_by must include mark")
        for column in hicrep_group_by:
            if column in text.columns and (text[column].str.strip() == "").any():
                errors.append(f"hicrep.group_by column {column!r} contains empty values")

    anchor_group_raw = config.get("anchor_consensus", {}).get("group_by", [])
    if not isinstance(anchor_group_raw, list) or not anchor_group_raw:
        errors.append(
            "anchor_consensus.group_by must be a non-empty list of assay-stratum columns"
        )
        anchor_group_by: list[str] = []
    else:
        anchor_group_by = [str(column) for column in anchor_group_raw]
        if len(anchor_group_by) != len(set(anchor_group_by)):
            errors.append("anchor_consensus.group_by must not contain duplicates")
        missing_anchor_columns = sorted(set(anchor_group_by) - set(text.columns))
        if missing_anchor_columns:
            errors.append(
                "anchor_consensus.group_by columns are absent from samples.tsv: "
                f"{missing_anchor_columns}"
            )
        required_anchor_columns = {
            "mark", "tissue", "library_protocol", "restriction_enzyme"
        }
        missing_required = sorted(required_anchor_columns - set(anchor_group_by))
        if missing_required:
            errors.append(
                "anchor_consensus.group_by must separate mark/tissue/protocol/enzyme "
                f"assay strata; missing: {missing_required}"
            )
        forbidden_anchor_columns = sorted(
            set(anchor_group_by)
            & {
                "sample_id", "subject_id", "replicate", "batch", "srr",
                "fastq_r1", "fastq_r2", "notes",
            }
        )
        if forbidden_anchor_columns:
            errors.append(
                "anchor_consensus.group_by must describe assay strata, not sample "
                f"identifiers or technical fields: {forbidden_anchor_columns}"
            )
        for column in anchor_group_by:
            if column in text.columns and (text[column].str.strip() == "").any():
                errors.append(
                    f"anchor_consensus.group_by column {column!r} contains empty values"
                )

        # A comparison must use one common anchor search space across both arms.
        # If a contrast-defining column entered this key, each arm could define a
        # different FitHiChIP hypothesis universe.
        for comparison in config.get("differential", {}).get("comparisons", []) or []:
            case_filter = comparison.get("case_filter", {}) or {}
            control_filter = comparison.get("control_filter", {}) or {}
            contrast_columns = {
                column
                for column in set(case_filter) | set(control_filter)
                if case_filter.get(column) != control_filter.get(column)
            }
            leaked = sorted(contrast_columns & set(anchor_group_by))
            if leaked:
                errors.append(
                    f"comparison {comparison.get('name')!r} changes {leaked}, but "
                    "those columns are in anchor_consensus.group_by; both arms "
                    "must share one anchor search space"
                )

    if text["sample_id"].duplicated().any():
        errors.append("sample_id values must be unique")
    biological_key = [
        "subject_id",
        *[column for column in hicrep_group_by if column in text.columns],
    ]
    duplicate_units = text.duplicated(biological_key, keep=False)
    if duplicate_units.any():
        units = text.loc[duplicate_units, biological_key].drop_duplicates().to_dict("records")
        errors.append(
            f"each {'/'.join(biological_key)} biological unit must occur once; "
            f"merge technical runs before analysis: {units}"
        )
    unsafe = [sid for sid in text["sample_id"] if not SAFE_ID.fullmatch(sid)]
    if unsafe:
        errors.append(f"sample_id values are unsafe for paths/wildcards: {unsafe}")
    demonstration_raw = config.get("reporting", {}).get(
        "demonstration_samples", []
    )
    if not isinstance(demonstration_raw, list):
        errors.append("reporting.demonstration_samples must be a list")
    else:
        demonstration = [str(value) for value in demonstration_raw]
        if len(demonstration) != len(set(demonstration)):
            errors.append("reporting.demonstration_samples must not contain duplicates")
        unknown_demonstration = sorted(set(demonstration) - set(text["sample_id"]))
        if unknown_demonstration:
            errors.append(
                "reporting.demonstration_samples contains unknown sample IDs: "
                f"{unknown_demonstration}"
            )
    unsafe_marks = sorted({mark for mark in text["mark"] if not SAFE_ID.fullmatch(mark)})
    if unsafe_marks:
        errors.append(f"mark values are unsafe for paths/wildcards: {unsafe_marks}")
    if "srr" in text.columns:
        accession_owners: dict[str, list[str]] = {}
        for row in text[["sample_id", "srr"]].itertuples(index=False):
            for accession in (part.strip() for part in row.srr.split(",")):
                if accession:
                    accession_owners.setdefault(accession, []).append(row.sample_id)
        invalid_accessions = sorted({
            accession
            for value in text["srr"]
            for accession in (part.strip() for part in value.split(","))
            if accession and not SRA_ACCESSION.fullmatch(accession)
        })
        if invalid_accessions:
            errors.append(
                "srr must contain only comma-separated SRR/ERR/DRR accessions: "
                f"{invalid_accessions}"
            )
        duplicate_accessions = {
            accession: owners
            for accession, owners in sorted(accession_owners.items())
            if len(owners) > 1
        }
        if duplicate_accessions:
            errors.append(
                "SRA accessions must be unique across and within sample rows; "
                f"duplicates: {duplicate_accessions}"
            )
    critical = [
        "subject_id", "cell_type", "tissue", "replicate", "mark", "fastq_r1",
        "fastq_r2", "batch", "library_protocol", "restriction_enzyme",
    ]
    for column in critical:
        if (text[column].str.strip() == "").any():
            errors.append(f"{column} must be present for every sample")
    for mate in ("fastq_r1", "fastq_r2"):
        unsafe_paths = [
            value for value in text[mate]
            if not SAFE_PATH.fullmatch(value) or ".." in Path(value).parts
        ]
        if unsafe_paths:
            errors.append(f"{mate} contains unsafe paths: {unsafe_paths}")
    same_mates = [
        str(row.sample_id)
        for row in text.itertuples(index=False)
        if Path(str(row.fastq_r1)) == Path(str(row.fastq_r2))
    ]
    if same_mates:
        errors.append(
            "fastq_r1 and fastq_r2 must be different files for every paired-end "
            f"library; identical mates found for: {same_mates}"
        )
    bad_replicates = pd.to_numeric(text["replicate"], errors="coerce").isna()
    if bad_replicates.any():
        errors.append("replicate must be numeric for every sample")

    assembly = config.get("assembly")
    if assembly not in genome_cfg:
        errors.append(f"assembly {assembly!r} is absent from genome.yaml")
    else:
        required_genome = {
            "fasta", "chromsizes", "gtf", "blacklist", "digest_bed", "enzyme",
            "macs3_genome_size",
        }
        absent = sorted(required_genome - set(genome_cfg[assembly]))
        if absent:
            errors.append(f"genome {assembly!r} lacks keys: {absent}")
        enzymes = set(text["restriction_enzyme"].str.strip())
        if enzymes != {str(genome_cfg[assembly].get("enzyme", ""))}:
            errors.append(
                f"sample restriction_enzyme values {sorted(enzymes)} do not match "
                f"the {assembly} digest enzyme {genome_cfg[assembly].get('enzyme')!r}"
            )

    fastq_dir = str(config.get("fastq_dir", "")).strip()
    if not fastq_dir:
        errors.append("fastq_dir must be configured")
    elif "srr" in text.columns:
        for row in text.itertuples(index=False):
            if not str(row.srr).strip():
                continue
            expected_r1 = Path(fastq_dir) / f"{row.sample_id}_R1.fastq.gz"
            expected_r2 = Path(fastq_dir) / f"{row.sample_id}_R2.fastq.gz"
            if Path(row.fastq_r1) != expected_r1 or Path(row.fastq_r2) != expected_r2:
                errors.append(
                    f"SRA sample {row.sample_id} FASTQ paths must be {expected_r1} and "
                    f"{expected_r2}, matching fastq_dir"
                )

    known_marks = set(config.get("macs3", {}).get("marks", {}))
    unknown_marks = sorted(set(text["mark"]) - known_marks)
    if unknown_marks:
        errors.append(f"no MACS3 mode is configured for marks: {unknown_marks}")

    supported_interactions = {"Peak-to-Peak", "Peak-to-NonPeak", "Peak-to-ALL", "ALL-to-ALL"}
    supported_backgrounds = {"Coverage_Bias"}
    if config.get("fithichip", {}).get("interaction_type") not in supported_interactions:
        errors.append("fithichip.interaction_type is not supported")
    if config.get("fithichip", {}).get("background_type") not in supported_backgrounds:
        errors.append(
            "fithichip.background_type must be Coverage_Bias; FitHiChIP ICE_Bias "
            "requires the separate HiC-Pro/iced executable contract, which this "
            "cooler-based workflow does not provide"
        )
    if config.get("differential", {}).get("method") != "pyDESeq2":
        errors.append("differential.method must be pyDESeq2")

    configured_kb = [int(k) for k in config["cooler"]["bin_sizes_kb"]]
    if not configured_kb or any(value <= 0 for value in configured_kb):
        errors.append("cooler.bin_sizes_kb must contain positive integers")
        resolutions: set[int] = set()
        base_resolution = 0
    else:
        resolutions = {value * 1000 for value in configured_kb}
        base_resolution = min(resolutions)
        if len(resolutions) != len(configured_kb):
            errors.append("cooler.bin_sizes_kb must not contain duplicates")
        incompatible = sorted(value for value in resolutions if value % base_resolution)
        if incompatible:
            errors.append(
                "every cooler resolution must be an integer multiple of the base "
                f"resolution {base_resolution}: {incompatible}"
            )
        if int(config["fithichip"]["bin_size"]) != base_resolution:
            errors.append(
                "fithichip.bin_size must equal the finest cooler resolution because "
                "FitHiChIP consumes the base .cool directly"
            )
    oracle_raw = config["cooler"].get("oracle_bin_sizes_kb", [])
    oracle_configured_kb: list[int] = []
    if not isinstance(oracle_raw, list) or not oracle_raw:
        errors.append("cooler.oracle_bin_sizes_kb must be a non-empty list")
    else:
        try:
            oracle_configured_kb = [int(value) for value in oracle_raw]
            if any(
                isinstance(value, bool)
                or float(value) != int(value)
                or int(value) <= 0
                for value in oracle_raw
            ):
                raise ValueError
        except (TypeError, ValueError):
            errors.append(
                "cooler.oracle_bin_sizes_kb must contain positive integers"
            )
            oracle_configured_kb = []
        if len(oracle_configured_kb) != len(set(oracle_configured_kb)):
            errors.append("cooler.oracle_bin_sizes_kb must not contain duplicates")

    required_resolutions = {
        int(config["fithichip"]["bin_size"]), int(config["apa"]["bin_size"]),
        int(config["hicrep"]["bin_size"]), 5_000, 10_000, 25_000, 100_000,
    }
    if config.get("stripes", {}).get("enabled"):
        required_resolutions.add(int(config["stripes"]["resolution"]))
    if config.get("mustache", {}).get("enabled"):
        required_resolutions.add(int(config["mustache"]["resolution"]))
    required_resolutions |= {value * 1000 for value in oracle_configured_kb}
    absent_resolutions = sorted(required_resolutions - resolutions)
    if absent_resolutions:
        errors.append(f"cooler.bin_sizes_kb lacks required bp resolutions: {absent_resolutions}")
    oracle_resolutions = sorted(value * 1000 for value in oracle_configured_kb)
    if any(
        coarse % fine
        for fine, coarse in zip(oracle_resolutions[:-1], oracle_resolutions[1:])
    ):
        errors.append(
            "successive cooler.oracle_bin_sizes_kb values must divide exactly so "
            "every fine ORACLE node has one deterministic coarse parent"
        )

    balance = config["cooler"].get("balance", {})
    if float(balance.get("tolerance", 0)) <= 0:
        errors.append("cooler.balance.tolerance must be positive")
    if int(balance.get("max_iterations", 0)) < 1:
        errors.append("cooler.balance.max_iterations must be at least one")
    if int(balance.get("min_nnz", 0)) < 0:
        errors.append("cooler.balance.min_nnz must be non-negative")
    if int(balance.get("ignore_diags", -1)) < 0:
        errors.append("cooler.balance.ignore_diags must be non-negative")

    fithichip = config["fithichip"]
    if not (0 < float(fithichip["fdr_threshold"]) <= 1):
        errors.append("fithichip.fdr_threshold must be in (0, 1]")
    if not (0 <= int(fithichip["lower_distance"]) < int(fithichip["upper_distance"])):
        errors.append("fithichip distances must satisfy 0 <= lower < upper")
    if int(fithichip.get("min_reads", 0)) < 0:
        errors.append("fithichip.min_reads must be non-negative")
    use_p2p = fithichip.get("use_p2p_background")
    if isinstance(use_p2p, bool) or not isinstance(use_p2p, int) or use_p2p not in {0, 1}:
        errors.append("fithichip.use_p2p_background must be the integer 0 or 1")
    if not isinstance(fithichip.get("merge_nearby"), bool):
        errors.append("fithichip.merge_nearby must be boolean")

    apa = config["apa"]
    apa_window = int(apa["window_size"])
    apa_bin = int(apa["bin_size"])
    if apa_window < 1 or apa_bin < 1:
        errors.append("apa.window_size and apa.bin_size must be positive")
    elif int(apa["min_loop_dist"]) < (2 * apa_window + 1) * apa_bin:
        errors.append(
            "apa.min_loop_dist must keep the aggregate window off the main diagonal"
        )
    if int(apa.get("n_random_controls", 0)) < 1:
        errors.append("apa.n_random_controls must be at least one")
    if float(apa.get("control_marginal_log2_tolerance", 0)) <= 0:
        errors.append("apa.control_marginal_log2_tolerance must be positive")
    if int(apa.get("max_control_attempts_per_draw", 0)) < 1:
        errors.append("apa.max_control_attempts_per_draw must be at least one")
    if int(apa.get("min_loops_for_apa", 0)) < 2:
        errors.append("apa.min_loops_for_apa must be at least two")
    apa_tolerance = apa.get("candidate_tolerance_bins")
    if (
        isinstance(apa_tolerance, bool)
        or not isinstance(apa_tolerance, int)
        or apa_tolerance != 1
    ):
        errors.append(
            "apa.candidate_tolerance_bins must be integer 1 to reconcile "
            "q-filtered donor calls without changing exact differential pixels"
        )

    differential = config.get("differential", {})
    if differential.get("hypothesis_source") != "fithichip_all_interactions":
        errors.append(
            "differential.hypothesis_source must be fithichip_all_interactions; "
            "q-filtered or merged calls outcome-select the tested universe"
        )
    if "min_sample_support" in differential:
        errors.append(
            "differential.min_sample_support is obsolete; use the explicit "
            "min_count and min_samples abundance contract"
        )
    tolerance = differential.get("candidate_tolerance_bins")
    if isinstance(tolerance, bool) or not isinstance(tolerance, int) or tolerance != 0:
        errors.append(
            "differential.candidate_tolerance_bins must be integer 0 because "
            "native FitHiChIP all-interaction pixels use one exact grid"
        )
    for key, floor in (("min_count", 1), ("min_samples", 1)):
        value = differential.get(key)
        if (
            isinstance(value, bool)
            or not isinstance(value, int)
            or value < floor
        ):
            errors.append(f"differential.{key} must be an integer >= {floor}")
    publication_min = differential.get("publication_min_complete_pairs")
    if (
        isinstance(publication_min, bool)
        or not isinstance(publication_min, int)
        or publication_min < 3
    ):
        errors.append(
            "differential.publication_min_complete_pairs must be an integer >= 3"
        )
        publication_min = 3
    require_publication = differential.get("require_publication_ready")
    if not isinstance(require_publication, bool):
        errors.append("differential.require_publication_ready must be boolean")
        require_publication = False
    if int(config.get("mustache", {}).get("comparison_tolerance_bins", -1)) < 0:
        errors.append("mustache.comparison_tolerance_bins must be non-negative")

    balance = config.get("cooler", {}).get("balance", {})
    if balance.get("weight_name", "weight") != "weight":
        errors.append(
            "cooler.balance.weight_name must be 'weight' because hicmatrix/"
            "pyGenomeTracks implicitly selects that column"
        )
    try:
        if float(balance.get("tolerance", 1e-5)) <= 0:
            errors.append("cooler.balance.tolerance must be positive")
    except (TypeError, ValueError):
        errors.append("cooler.balance.tolerance must be numeric")
    try:
        max_iterations = int(balance.get("max_iterations", 200))
        if max_iterations <= 0 or max_iterations != float(balance.get("max_iterations", 200)):
            errors.append("cooler.balance.max_iterations must be a positive integer")
    except (TypeError, ValueError):
        errors.append("cooler.balance.max_iterations must be a positive integer")

    paired_by = differential.get("paired_by")
    comparisons = differential.get("comparisons", []) or []
    if comparisons and (
        not isinstance(paired_by, str) or not paired_by.strip()
    ):
        errors.append(
            "differential.paired_by must name a non-empty sample-sheet column "
            "when comparisons are configured; this release supports paired "
            "differential inference only"
        )
        paired_by = None
    names = [c.get("name", "") for c in comparisons]
    if len(names) != len(set(names)) or any(not SAFE_ID.fullmatch(name) for name in names):
        errors.append("differential comparison names must be unique safe identifiers")
    for comp in comparisons:
        include_subjects_raw = comp.get("include_subjects")
        include_subjects: list[str] | None = None
        if include_subjects_raw is not None:
            if not isinstance(include_subjects_raw, list) or not include_subjects_raw:
                errors.append(
                    f"comparison {comp.get('name')} include_subjects must be a "
                    "non-empty list"
                )
            else:
                include_subjects = [str(value).strip() for value in include_subjects_raw]
                if any(not value for value in include_subjects):
                    errors.append(
                        f"comparison {comp.get('name')} include_subjects contains "
                        "an empty value"
                    )
                if len(include_subjects) != len(set(include_subjects)):
                    errors.append(
                        f"comparison {comp.get('name')} include_subjects contains duplicates"
                    )
                unknown = sorted(set(include_subjects) - set(text["subject_id"]))
                if unknown:
                    errors.append(
                        f"comparison {comp.get('name')} include_subjects contains "
                        f"unknown subject IDs: {unknown}"
                    )
        try:
            cases = _selected(
                text, comp.get("case_filter", {}), comp.get("mark"), include_subjects
            )
            controls = _selected(
                text, comp.get("control_filter", {}), comp.get("mark"), include_subjects
            )
        except ValueError as exc:
            errors.append(f"comparison {comp.get('name')}: {exc}")
            continue
        overlap = set(cases["sample_id"]) & set(controls["sample_id"])
        if overlap:
            errors.append(f"comparison {comp['name']} overlaps case/control samples: {sorted(overlap)}")
        if len(cases) < 2 or len(controls) < 2:
            errors.append(f"comparison {comp['name']} needs at least two biological samples per arm")
        combined = pd.concat([cases, controls])
        configured_min_samples = differential.get("min_samples")
        if (
            isinstance(configured_min_samples, int)
            and not isinstance(configured_min_samples, bool)
            and configured_min_samples > len(combined)
        ):
            errors.append(
                f"comparison {comp['name']} selects {len(combined)} samples, below "
                f"differential.min_samples={configured_min_samples}"
            )
        for column in ("mark", "tissue", "library_protocol"):
            if combined[column].nunique() > 1:
                errors.append(f"comparison {comp['name']} mixes {column} values")
        covariates = config.get("differential", {}).get("covariates", []) or []
        for covariate in covariates:
            if covariate not in text.columns:
                errors.append(f"comparison {comp['name']} covariate {covariate!r} is absent")
        if combined["batch"].nunique() > 1 and "batch" not in covariates:
            errors.append(
                f"comparison {comp['name']} mixes batch values but batch is not "
                "listed in differential.covariates"
            )
        if paired_by:
            if paired_by not in text.columns:
                errors.append(f"differential.paired_by {paired_by!r} is absent from samples.tsv")
            elif set(cases[paired_by]) != set(controls[paired_by]):
                errors.append(f"comparison {comp['name']} has unmatched {paired_by} levels")
            elif cases[paired_by].duplicated().any() or controls[paired_by].duplicated().any():
                errors.append(f"comparison {comp['name']} is not one-to-one paired by {paired_by}")
            elif require_publication and len(set(cases[paired_by])) < publication_min:
                errors.append(
                    f"comparison {comp['name']} has {len(set(cases[paired_by]))} "
                    f"complete {paired_by} pairs, below publication minimum "
                    f"{publication_min}; set require_publication_ready=false to "
                    "run it explicitly as PILOT_UNDERPOWERED"
                )

        # Reject exact confounding before pyDESeq2 starts. Strings are categorical,
        # matching the sample-sheet contract and pyDESeq2 design_factors behavior.
        design = combined.copy()
        design["condition"] = [
            "case" if sid in set(cases["sample_id"]) else "control"
            for sid in design["sample_id"]
        ]
        factors = [factor for factor in ([paired_by] if paired_by else []) + covariates
                   if factor in design.columns and design[factor].nunique() > 1]
        factors.append("condition")
        matrix = pd.get_dummies(
            design[list(dict.fromkeys(factors))].astype(str), drop_first=True, dtype=float
        )
        matrix.insert(0, "intercept", 1.0)
        if np.linalg.matrix_rank(matrix.to_numpy()) < matrix.shape[1]:
            errors.append(
                f"comparison {comp['name']} has a rank-deficient design; a covariate "
                "is confounded with condition or pairing"
            )

    viz = config.get("viz", {})
    if viz.get("assembly") != assembly:
        errors.append(
            f"viz.assembly {viz.get('assembly')!r} does not match assembly {assembly!r}; "
            "replace the example loci for the selected assembly"
        )
    regions = viz.get("regions", [])
    region_names = [str(region.get("name", "")) for region in regions]
    if len(region_names) != len(set(region_names)) or any(
        not SAFE_ID.fullmatch(name) for name in region_names
    ):
        errors.append("viz region names must be unique safe identifiers")
    for region in regions:
        try:
            start, end, viewpoint = int(region["start"]), int(region["end"]), int(region["viewpoint"])
            if not (0 <= start < viewpoint < end):
                errors.append(f"viz region {region.get('name')} viewpoint must lie inside [start,end)")
        except (KeyError, TypeError, ValueError):
            errors.append(f"viz region {region.get('name')} needs integer start/end/viewpoint")
        if not SAFE_ID.fullmatch(str(region.get("chrom", ""))):
            errors.append(f"viz region {region.get('name')} needs a safe chromosome name")
        if not str(region.get("viewpoint_label", "")).strip():
            errors.append(f"viz region {region.get('name')} needs viewpoint_label")

    oracle_export = config.get("oracle_export", {})
    if not isinstance(oracle_export.get("primary_chromosomes_only", True), bool):
        errors.append("oracle_export.primary_chromosomes_only must be boolean")
    drop_raw = oracle_export.get("drop_chromosomes", [])
    if not isinstance(drop_raw, list):
        errors.append("oracle_export.drop_chromosomes must be a list")
    else:
        drop_chromosomes = [str(value).strip() for value in drop_raw]
        if any(not value or not SAFE_ID.fullmatch(value) for value in drop_chromosomes):
            errors.append(
                "oracle_export.drop_chromosomes must contain non-empty safe "
                "chromosome names"
            )
        if len(drop_chromosomes) != len(set(drop_chromosomes)):
            errors.append("oracle_export.drop_chromosomes must not contain duplicates")
    if not isinstance(oracle_export.get("microbiome_metadata_tsv", ""), str):
        errors.append("oracle_export.microbiome_metadata_tsv must be a string path")

    if errors:
        raise ValueError("Configuration validation failed:\n- " + "\n- ".join(errors))
