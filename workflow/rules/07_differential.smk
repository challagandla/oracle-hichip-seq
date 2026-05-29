# Stage 07 — Differential loop analysis
# Compare tumor vs normal (or any two conditions) on the union loop set.
# Uses pyDESeq2 by default; switch to diffHic (R) via config["differential"]["method"].

def _group_samples_for_diff():
    """Define case/control groups from samples.tsv. Tumor vs healthy by default."""
    cases    = SAMPLES[SAMPLES["tissue"] == "tumor"]["sample_id"].tolist()
    controls = SAMPLES[SAMPLES["disease"] == "healthy"]["sample_id"].tolist()
    return cases, controls

CASES, CONTROLS = _group_samples_for_diff()

rule build_union_loops:
    """Merge all per-sample FitHiChIP loops into a union BEDPE."""
    input:
        loops = expand(RESULTS / "loops/{sample}/{sample}.interactions_FitHiC_Q0.01.bed",
                       sample=CASES + CONTROLS)
    output:
        bedpe = RESULTS / "diff/union_loops.bedpe"
    log:
        RESULTS / "logs/build_union/union.log"
    shell:
        r"""
        cat {input.loops} | \
            awk 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,$4,$5,$6}}' | \
            sort -k1,1 -k2,2n -k4,4 -k5,5n | uniq > {output.bedpe} 2> {log}
        """


rule count_per_loop:
    """
    For each sample, count valid pairs supporting each union-loop. Implements
    a fast cooler-based pile-up at loop coordinates.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        bedpe = RESULTS / "diff/union_loops.bedpe"
    output:
        counts = RESULTS / "diff/counts/{sample}.counts.tsv"
    params:
        res = config["fithichip"]["bin_size"]
    threads: 4
    log:
        RESULTS / "logs/count_per_loop/{sample}.log"
    script:
        "../scripts/count_per_loop.py"


rule differential_loops:
    """Run pyDESeq2 (or diffHic) and emit a results TSV + volcano PNG."""
    input:
        counts = expand(RESULTS / "diff/counts/{sample}.counts.tsv",
                        sample=CASES + CONTROLS)
    output:
        tsv = RESULTS / "diff/differential_loops.tsv",
        volcano = RESULTS / "diff/volcano.png",
        ma = RESULTS / "diff/ma_plot.png"
    params:
        cases = CASES,
        controls = CONTROLS,
        method = config["differential"]["method"],
        fdr = config["differential"]["fdr"],
        log2fc_min = config["differential"]["log2fc_min"]
    threads: 8
    log:
        RESULTS / "logs/differential_loops/diff.log"
    script:
        "../scripts/differential_loops.py"
