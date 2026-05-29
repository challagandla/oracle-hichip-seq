# Stage 06 — Matrix and loop QC
# Cis/trans ratio, P(s) distance decay, insulation scores, A/B compartments,
# APA aggregate analysis on the called loop set, HiCRep replicate concordance.

rule cooltools_expected_cis:
    """
    P(s) distance-decay curve. The expected −1 slope on log–log is the
    canonical sanity check for any HiC-style assay.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool"
    output:
        expected = RESULTS / "qc/expected/{sample}.expected.cis.tsv"
    params:
        res = 25000
    threads: 4
    log:
        RESULTS / "logs/cooltools_expected/{sample}.log"
    shell:
        r"""
        cooltools expected-cis -p {threads} \
            {input.mcool}::resolutions/{params.res} \
            > {output.expected} 2> {log}
        """


rule cooltools_insulation:
    """Insulation score (TAD boundary detection)."""
    input:
        mcool = RESULTS / "cool/{sample}.mcool"
    output:
        tsv = RESULTS / "qc/insulation/{sample}.insulation.tsv"
    params:
        res = 25000,
        window = 250000
    threads: 4
    log:
        RESULTS / "logs/cooltools_insulation/{sample}.log"
    shell:
        r"""
        cooltools insulation \
            -p {threads} \
            {input.mcool}::resolutions/{params.res} \
            {params.window} \
            > {output.tsv} 2> {log}
        """


rule cooltools_eigs_cis:
    """A/B compartment eigenvectors at 100 kb."""
    input:
        mcool = RESULTS / "cool/{sample}.mcool"
    output:
        cis = RESULTS / "qc/compartments/{sample}.cis.eigs.tsv"
    params:
        res = 100000
    threads: 4
    log:
        RESULTS / "logs/cooltools_eigs/{sample}.log"
    shell:
        r"""
        cooltools eigs-cis -p {threads} \
            {input.mcool}::resolutions/{params.res} \
            -o {output.cis} 2> {log}
        """


rule compartments_to_bigwig:
    """Export E1 (A/B compartments) to a bigWig for browser / pyGenomeTracks."""
    input:
        eigs = RESULTS / "qc/compartments/{sample}.cis.eigs.tsv"
    output:
        bw = RESULTS / "qc/compartments/{sample}.E1.bw"
    params:
        chromsizes = GENOME["chromsizes"]
    threads: 1
    log:
        RESULTS / "logs/compartments_to_bigwig/{sample}.log"
    script:
        "../scripts/compartments_to_bigwig.py"


rule hicrep_replicate_qc:
    """
    Stratum-adjusted correlation between biological replicates of the same
    subject + mark. Only metric robust to distance-decay differences.
    """
    input:
        mcools = lambda wc: expand(
            RESULTS / "cool/{sample}.mcool",
            sample=SAMPLES[
                (SAMPLES["subject_id"] == SAMPLES.loc[wc.sample, "subject_id"]) &
                (SAMPLES["mark"] == SAMPLES.loc[wc.sample, "mark"])
            ]["sample_id"].tolist()
        )
    output:
        json = RESULTS / "qc/hicrep/{sample}.hicrep.json"
    params:
        bin = config["hicrep"]["bin_size"],
        maxd = config["hicrep"]["max_dist"],
        h = config["hicrep"]["h_smooth"]
    threads: 4
    log:
        RESULTS / "logs/hicrep/{sample}.log"
    script:
        "../scripts/hicrep_replicate_qc.py"


rule apa_plot:
    """
    Aggregate Peak Analysis on the called loop set with random-shift controls.
    A score ≥ 1.5 vs controls is the standard quality cutoff.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        loops = RESULTS / "loops/{sample}/{sample}.interactions_FitHiC_Q0.01.bed"
    output:
        png = RESULTS / "qc/apa/{sample}.apa.png",
        json = RESULTS / "qc/apa/{sample}.apa.json"
    params:
        window = config["apa"]["window_size"],
        bin_size = config["apa"]["bin_size"],
        min_dist = config["apa"]["min_loop_dist"],
        n_ctrl = config["apa"]["n_random_controls"]
    threads: 4
    log:
        RESULTS / "logs/apa/{sample}.log"
    script:
        "../scripts/apa_plot.py"


rule loop_qc_summary:
    """
    Aggregate every QC metric into a single JSON per sample with pass/fail flags.
    Consumed by MultiQC custom content.
    """
    input:
        pair_stats = RESULTS / "qc/pairtools/{sample}.pairs.stats.txt",
        dedup_stats = RESULTS / "qc/pairtools/{sample}.dedup.stats.txt",
        expected = RESULTS / "qc/expected/{sample}.expected.cis.tsv",
        apa_json = RESULTS / "qc/apa/{sample}.apa.json",
        hicrep = RESULTS / "qc/hicrep/{sample}.hicrep.json",
        loops = RESULTS / "loops/{sample}/{sample}.interactions_FitHiC_Q0.01.bed"
    output:
        json = RESULTS / "qc/loop_qc/{sample}.json",
        md   = RESULTS / "qc/loop_qc/{sample}.md"
    threads: 1
    log:
        RESULTS / "logs/loop_qc_summary/{sample}.log"
    script:
        "../scripts/loop_qc_summary.py"
