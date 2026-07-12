# Stage 11 — Cohort-level publication figures
#
# The per-sample PNGs the earlier stages emit are diagnostics: one APA per library,
# one browser panel per locus. Nothing assembled them across the cohort, so the
# questions a reader actually asks -- is the loop count driven by biology or by
# depth, do the subsets separate, are stripes an artefact of the anchor mark --
# had no figure that answered them. This stage produces those.
#
# Figures depend on the QC/loop/diff artefacts, not on the raw data, so they are
# cheap to re-render after a threshold change.


rule publication_figures:
    input:
        loops   = expand(RESULTS / f"loops/{{sample}}/{{sample}}.interactions_FitHiC_{FITHICHIP_Q_LABEL}.bed",
                         sample=SAMPLE_IDS),
        apa     = expand(RESULTS / "qc/apa/{sample}.apa.npy", sample=SAMPLE_IDS),
        hicrep  = expand(RESULTS / "qc/hicrep/{sample}.hicrep.json", sample=SAMPLE_IDS),
        eigs    = expand(RESULTS / "qc/compartments/{sample}.cis.eigs.tsv", sample=SAMPLE_IDS),
        expect  = expand(RESULTS / "qc/expected/{sample}.expected.cis.tsv", sample=SAMPLE_IDS),
        dedup   = expand(RESULTS / "qc/pairtools/{sample}.dedup.stats.txt", sample=SAMPLE_IDS),
        stripes = (expand(RESULTS / "stripes/{sample}/result_filtered.tsv", sample=SAMPLE_IDS)
                   if config.get("stripes", {}).get("enabled") else []),
        diff    = expand(RESULTS / "diff/{comparison}/differential_loops.tsv",
                         comparison=[c["name"] for c in config.get("differential", {}).get("comparisons", [])]),
    output:
        table = RESULTS / "figures/library_summary.tsv",
        figs = expand(RESULTS / "figures/{fig}.{ext}",
                      fig=["figure1_library_qc", "figure2_reproducibility",
                           "figure3_loops_apa", "figure4_differential",
                           "figure5_stripes"],
                      ext=["pdf", "png"]),
    params:
        results     = lambda wc: str(RESULTS),
        outdir      = lambda wc: str(RESULTS / "figures"),
        samples_tsv = config["samples_tsv"],
        q_label     = FITHICHIP_Q_LABEL,
    threads: 2
    conda: "../envs/figures.yaml"
    log:
        RESULTS / "logs/figures/publication_figures.log"
    script:
        "../scripts/figures.py"
