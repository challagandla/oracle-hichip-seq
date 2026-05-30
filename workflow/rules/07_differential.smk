# Stage 07 — Differential loop analysis
# Compare case/control groups on a union loop set. Comparisons are explicit in
# config.yaml and must be mark/tissue/protocol compatible.

import json

DIFF_COMPARISONS = config.get("differential", {}).get("comparisons", []) or []


def _comparison_cfg(name):
    for comp in DIFF_COMPARISONS:
        if comp["name"] == name:
            return comp
    raise ValueError(f"Differential comparison {name!r} not found in config['differential']['comparisons']")


def _filter_samples(filter_dict, mark=None):
    df = SAMPLES.copy()
    if mark is not None:
        df = df[df["mark"] == mark]
    for key, value in (filter_dict or {}).items():
        if key not in df.columns:
            raise ValueError(f"Differential filter column {key!r} is not in samples.tsv")
        df = df[df[key] == str(value)]
    return df["sample_id"].tolist()


def _group_samples_for_diff(wc):
    comp = _comparison_cfg(wc.comparison)
    mark = comp.get("mark")
    cases = _filter_samples(comp.get("case_filter", {}), mark=mark)
    controls = _filter_samples(comp.get("control_filter", {}), mark=mark)
    return cases, controls, comp


def _loop_inputs_for_comparison(wc):
    cases, controls, comp = _group_samples_for_diff(wc)
    if len(cases) == 0 or len(controls) == 0:
        raise ValueError(
            f"Comparison {wc.comparison!r} has cases={cases} and controls={controls}. "
            "Define matched, same-mark case/control samples in config.yaml."
        )
    return expand(
        RESULTS / f"loops/{{sample}}/{{sample}}.interactions_FitHiC_{FITHICHIP_Q_LABEL}.bed",
        sample=cases + controls,
    )


def _count_inputs_for_comparison(wc):
    cases, controls, comp = _group_samples_for_diff(wc)
    return expand(RESULTS / "diff/{comparison}/counts/{sample}.counts.tsv",
                  comparison=wc.comparison, sample=cases + controls)


rule build_union_loops:
    """Merge all per-sample FitHiChIP loops into a comparison-specific union BEDPE."""
    input:
        loops = _loop_inputs_for_comparison
    output:
        bedpe = RESULTS / "diff/{comparison}/union_loops.bedpe"
    log:
        RESULTS / "logs/build_union/{comparison}.log"
    shell:
        r"""
        cat {input.loops} | \
            awk 'BEGIN{{OFS="\t"}} !/^#/ {{print $1,$2,$3,$4,$5,$6}}' | \
            sort -k1,1 -k2,2n -k4,4 -k5,5n | uniq > {output.bedpe} 2> {log}
        test -s {output.bedpe}
        """


rule count_per_loop:
    """
    For each sample, count valid pairs supporting each union-loop using the
    unbalanced cooler matrix at the FitHiChIP bin size.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        bedpe = RESULTS / "diff/{comparison}/union_loops.bedpe"
    output:
        counts = RESULTS / "diff/{comparison}/counts/{sample}.counts.tsv"
    params:
        res = config["fithichip"]["bin_size"]
    threads: 4
    log:
        RESULTS / "logs/count_per_loop/{comparison}_{sample}.log"
    script:
        "../scripts/count_per_loop.py"


rule differential_loops:
    """Run pyDESeq2 (or diffHic later) and emit a results TSV + volcano PNG."""
    input:
        counts = _count_inputs_for_comparison
    output:
        tsv = RESULTS / "diff/{comparison}/differential_loops.tsv",
        volcano = RESULTS / "diff/{comparison}/volcano.png",
        ma = RESULTS / "diff/{comparison}/ma_plot.png",
        design = RESULTS / "diff/{comparison}/design.json"
    params:
        groups = lambda wc: _group_samples_for_diff(wc),
        method = config["differential"]["method"],
        fdr = config["differential"]["fdr"],
        log2fc_min = config["differential"]["log2fc_min"]
    threads: 8
    log:
        RESULTS / "logs/differential_loops/{comparison}.log"
    script:
        "../scripts/differential_loops.py"
