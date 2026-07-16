# Final, machine-readable provenance after all reportable analyses complete.

PROVENANCE_REFERENCE_CONTRACT = {
    "fasta": GENOME["fasta"],
    "chromsizes": GENOME["chromsizes"],
    "gtf": GENOME["gtf"],
    "blacklist": GENOME["blacklist"],
    "restriction_digest_bed": GENOME["digest_bed"],
}
PROVENANCE_REFERENCE_CONTRACT.update({
    f"alignment_index_{i + 1}_{Path(path).name}": path
    for i, path in enumerate(BWA_INDEX_FILES)
})

rule provenance_manifest:
    input:
        samples = config["samples_tsv"],
        parameters = "config/config.yaml",
        genome = config["genome_yaml"],
        runner = "environment.runner.yml",
        envs = sorted(str(path) for path in Path("workflow/envs").glob("*.yaml")),
        reference_fasta = GENOME["fasta"],
        reference_chromsizes = GENOME["chromsizes"],
        reference_gtf = GENOME["gtf"],
        reference_blacklist = GENOME["blacklist"],
        reference_digest = GENOME["digest_bed"],
        reference_alignment_indexes = BWA_INDEX_FILES,
        multiqc = RESULTS / "multiqc/multiqc_report.html",
        figures = RESULTS / "figures/library_summary.tsv",
        figure_files = expand(
            RESULTS / "figures/{fig}.{ext}",
            fig=["figure1_library_qc", "figure2_reproducibility",
                 "figure3_loops_apa", "figure4_differential", "figure5_stripes"],
            ext=["pdf", "png"],
        ),
        loop_qc = expand(RESULTS / "qc/loop_qc/{sample}.json", sample=SAMPLE_IDS),
        balance_qc = expand(RESULTS / "qc/balance/{sample}.balance.json", sample=SAMPLE_IDS),
        blacklist_qc = expand(RESULTS / "qc/blacklist/{sample}.blacklist_filter.json", sample=SAMPLE_IDS),
        loop_call_audits = expand(
            RESULTS / "loops/{sample}/{sample}.interactions_FitHiC.all.audit.json",
            sample=SAMPLE_IDS,
        ),
        oracle = expand(RESULTS / "oracle_cos/{sample}.manifest.json", sample=SAMPLE_IDS),
        stripes = ([RESULTS / "stripes/stripe_summary.tsv"]
                   if config.get("stripes", {}).get("enabled") else []),
        differential = (
            expand(
                RESULTS / "diff/{comparison}/design.json",
                comparison=[c["name"] for c in config.get("differential", {}).get("comparisons", [])],
            )
            + expand(
                RESULTS / "diff/{comparison}/paired_effects.tsv",
                comparison=[c["name"] for c in config.get("differential", {}).get("comparisons", [])],
            )
        ),
        hypothesis_universes = expand(
            RESULTS / "diff/{comparison}/hypothesis_universe.json",
            comparison=[c["name"] for c in config.get("differential", {}).get("comparisons", [])],
        ),
    output:
        json = RESULTS / "provenance/run_manifest.json",
    params:
        repository = ".",
        effective_config = config,
        reference_contract = PROVENANCE_REFERENCE_CONTRACT,
        conda_cache = ".snakemake/conda",
        runner_environment = "oracle-hichip-runner",
    conda: "../envs/pandas.yaml"
    log:
        RESULTS / "logs/provenance.log",
    script:
        "../scripts/provenance_manifest.py"
