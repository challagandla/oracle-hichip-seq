# Stage 07 — Differential loop analysis
# Compare case/control groups on an exact-grid, q-unthresholded FitHiChIP
# hypothesis universe. Comparisons are explicit in config.yaml and must be
# mark/tissue/protocol compatible.

DIFF_COMPARISONS = config.get("differential", {}).get("comparisons", []) or []


def _comparison_cfg(name):
    for comp in DIFF_COMPARISONS:
        if comp["name"] == name:
            return comp
    raise ValueError(f"Differential comparison {name!r} not found in config['differential']['comparisons']")


def _filter_samples(filter_dict, mark=None, include_subjects=None):
    df = SAMPLES.copy()
    if mark is not None:
        df = df[df["mark"] == mark]
    for key, value in (filter_dict or {}).items():
        if key not in df.columns:
            raise ValueError(f"Differential filter column {key!r} is not in samples.tsv")
        df = df[df[key] == str(value)]
    if include_subjects is not None:
        df = df[df["subject_id"].isin([str(value) for value in include_subjects])]
    return df["sample_id"].tolist()


def _group_samples_for_diff(wc):
    comp = _comparison_cfg(wc.comparison)
    mark = comp.get("mark")
    include_subjects = comp.get("include_subjects")
    cases = _filter_samples(
        comp.get("case_filter", {}), mark=mark,
        include_subjects=include_subjects,
    )
    controls = _filter_samples(
        comp.get("control_filter", {}), mark=mark,
        include_subjects=include_subjects,
    )
    return cases, controls, comp


def _samples_for_comparison(wc):
    cases, controls, comp = _group_samples_for_diff(wc)
    if len(cases) == 0 or len(controls) == 0:
        raise ValueError(
            f"Comparison {wc.comparison!r} has cases={cases} and controls={controls}. "
            "Define matched, same-mark case/control samples in config.yaml."
        )
    return cases + controls


def _all_interaction_inputs_for_comparison(wc):
    return expand(
        RESULTS / "loops/{sample}/{sample}.interactions_FitHiC.all.tsv.gz",
        sample=_samples_for_comparison(wc),
    )


def _eligible_inputs_for_comparison(wc):
    return expand(
        RESULTS / "loops/{sample}/{sample}.interactions_FitHiC.eligible.tsv.gz",
        sample=_samples_for_comparison(wc),
    )


def _normalization_audits_for_comparison(wc):
    return expand(
        RESULTS / "loops/{sample}/{sample}.interactions_FitHiC.all.audit.json",
        sample=_samples_for_comparison(wc),
    )


def _count_inputs_for_comparison(wc):
    cases, controls, comp = _group_samples_for_diff(wc)
    return expand(RESULTS / "diff/{comparison}/counts/{sample}.counts.tsv",
                  comparison=wc.comparison, sample=cases + controls)


rule build_union_loops:
    """Create a label-blind exact-grid universe from unthresholded interactions."""
    input:
        eligible = _eligible_inputs_for_comparison,
        all_interactions = _all_interaction_inputs_for_comparison,
        audits = _normalization_audits_for_comparison,
        blacklist = GENOME["blacklist"],
        chromsizes = GENOME["chromsizes"],
        normalizer_code = "workflow/scripts/normalize_fithichip_all.py",
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        bedpe = RESULTS / "diff/{comparison}/union_loops.bedpe",
        support = RESULTS / "diff/{comparison}/candidate_support.tsv",
        manifest = RESULTS / "diff/{comparison}/hypothesis_universe.json",
    params:
        bin_size = config["fithichip"]["bin_size"],
        lower_distance = config["fithichip"]["lower_distance"],
        upper_distance = config["fithichip"]["upper_distance"],
        min_count = config["differential"].get("min_count", 5),
        min_samples = config["differential"].get("min_samples", 2),
        samples = lambda wc: _samples_for_comparison(wc),
    conda: "../envs/pandas.yaml"
    log:
        RESULTS / "logs/build_union/{comparison}.log"
    script:
        "../scripts/build_differential_universe.py"


rule count_per_loop:
    """
    For each sample, count valid pairs supporting each union-loop using the
    unbalanced cooler matrix at the FitHiChIP bin size.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        bedpe = RESULTS / "diff/{comparison}/union_loops.bedpe",
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        counts = RESULTS / "diff/{comparison}/counts/{sample}.counts.tsv"
    params:
        res = config["fithichip"]["bin_size"],
        lower_distance = config["fithichip"]["lower_distance"],
        upper_distance = config["fithichip"]["upper_distance"],
    threads: 4
    conda: "../envs/coolerpy.yaml"
    log:
        RESULTS / "logs/count_per_loop/{comparison}_{sample}.log"
    script:
        "../scripts/count_per_loop.py"


rule differential_loops:
    """Run the validated paired pyDESeq2 model and emit tables and QC plots."""
    input:
        counts = _count_inputs_for_comparison,
        universe = RESULTS / "diff/{comparison}/union_loops.bedpe",
        support = RESULTS / "diff/{comparison}/candidate_support.tsv",
        universe_manifest = RESULTS / "diff/{comparison}/hypothesis_universe.json",
        source_audits = _normalization_audits_for_comparison,
        samples = config["samples_tsv"],
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        tsv = RESULTS / "diff/{comparison}/differential_loops.tsv",
        volcano = RESULTS / "diff/{comparison}/volcano.png",
        ma = RESULTS / "diff/{comparison}/ma_plot.png",
        design = RESULTS / "diff/{comparison}/design.json",
        paired_effects = RESULTS / "diff/{comparison}/paired_effects.tsv",
    params:
        groups = lambda wc: _group_samples_for_diff(wc),
        method = config["differential"]["method"],
        fdr = config["differential"]["fdr"],
        log2fc_min = config["differential"]["log2fc_min"],
        paired_by = config["differential"].get("paired_by", "subject_id"),
        covariates = config["differential"].get("covariates", []),
        min_count = config["differential"].get("min_count", 5),
        min_samples = config["differential"].get("min_samples", 2),
        publication_min_complete_pairs = config["differential"].get(
            "publication_min_complete_pairs", 3
        ),
        require_publication_ready = config["differential"].get(
            "require_publication_ready", False
        ),
    threads: 8
    conda: "../envs/pydeseq2.yaml"
    log:
        RESULTS / "logs/differential_loops/{comparison}.log"
    script:
        "../scripts/differential_loops.py"
