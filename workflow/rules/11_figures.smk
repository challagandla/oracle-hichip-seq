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
        samples = config["samples_tsv"],
        loops   = expand(RESULTS / f"loops/{{sample}}/{{sample}}.interactions_FitHiC_{FITHICHIP_Q_LABEL}.bed",
                         sample=SAMPLE_IDS),
        apa     = expand(RESULTS / "qc/apa/{sample}.apa.npy", sample=SAMPLE_IDS),
        apa_json = expand(RESULTS / "qc/apa/{sample}.apa.json", sample=SAMPLE_IDS),
        balance = expand(RESULTS / "qc/balance/{sample}.balance.json", sample=SAMPLE_IDS),
        loop_qc = expand(RESULTS / "qc/loop_qc/{sample}.json", sample=SAMPLE_IDS),
        hicrep  = expand(RESULTS / "qc/hicrep/{sample}.hicrep.json", sample=SAMPLE_IDS),
        contact_depth = expand(RESULTS / "qc/contact_depth/{sample}.json", sample=SAMPLE_IDS),
        eigs    = expand(RESULTS / "qc/compartments/{sample}.cis.eigs.tsv", sample=SAMPLE_IDS),
        expect  = expand(RESULTS / "qc/expected/{sample}.expected.cis.tsv", sample=SAMPLE_IDS),
        dedup   = expand(RESULTS / "qc/pairtools/{sample}.dedup.stats.txt", sample=SAMPLE_IDS),
        pairs   = expand(RESULTS / "qc/pairtools/{sample}.pairs.stats.txt", sample=SAMPLE_IDS),
        stripes = (expand(RESULTS / "stripes/{sample}/result_filtered.tsv", sample=SAMPLE_IDS)
                   if config.get("stripes", {}).get("enabled") else []),
        diff    = expand(RESULTS / "diff/{comparison}/differential_loops.tsv",
                         comparison=[c["name"] for c in config.get("differential", {}).get("comparisons", [])]),
        diff_design = expand(
            RESULTS / "diff/{comparison}/design.json",
            comparison=[c["name"] for c in config.get("differential", {}).get("comparisons", [])],
        ),
        paired_effects = expand(
            RESULTS / "diff/{comparison}/paired_effects.tsv",
            comparison=[c["name"] for c in config.get("differential", {}).get("comparisons", [])],
        ),
        shared_code = SHARED_SCRIPT_DEPS,
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
        q_label     = FITHICHIP_Q_LABEL,
        min_contacts = config["hicrep"]["min_contacts_for_scc"],
        hicrep_threshold = config["hicrep"]["threshold_pass"],
        differential_fdr = config["differential"]["fdr"],
        differential_log2fc_min = config["differential"]["log2fc_min"],
        qc_thresholds = config["qc_thresholds"],
        apa_bin_size = config["apa"]["bin_size"],
        comparisons = [
            c["name"] for c in config.get("differential", {}).get("comparisons", [])
        ],
        demonstration_samples = config.get("reporting", {}).get(
            "demonstration_samples", []
        ),
    threads: 2
    conda: "../envs/figures.yaml"
    log:
        RESULTS / "logs/figures/publication_figures.log"
    script:
        "../scripts/figures.py"
