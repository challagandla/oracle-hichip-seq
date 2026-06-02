# Stage 10 — Aggregate report
# MultiQC pulls FastQC, fastp, pairtools, MACS2, cooler, FitHiChIP, HiCRep,
# APA, and ORACLE QC into a single HTML.

rule multiqc_report:
    input:
        fastqc = expand(RESULTS / "qc/fastqc_raw/{sample}_R1_fastqc.html", sample=SAMPLE_IDS),
        fastp  = expand(RESULTS / "qc/fastp/{sample}.fastp.json", sample=SAMPLE_IDS),
        pairs  = expand(RESULTS / "qc/pairtools/{sample}.pairs.stats.txt", sample=SAMPLE_IDS),
        dedup  = expand(RESULTS / "qc/pairtools/{sample}.dedup.stats.txt", sample=SAMPLE_IDS),
        apa    = expand(RESULTS / "qc/apa/{sample}.apa.json", sample=SAMPLE_IDS),
        loop   = expand(RESULTS / "qc/loop_qc/{sample}.json", sample=SAMPLE_IDS),
        hicrep = expand(RESULTS / "qc/hicrep/{sample}.hicrep.json", sample=SAMPLE_IDS)
    output:
        html = RESULTS / "multiqc/multiqc_report.html"
    params:
        outdir    = RESULTS / "multiqc",
        searchdir = RESULTS,
        mqc_cfg   = "config/multiqc_config.yaml"
    conda: "../envs/multiqc.yaml"
    log:
        RESULTS / "logs/multiqc/multiqc.log"
    shell:
        r"""
        multiqc -f \
            --config {params.mqc_cfg} \
            -o {params.outdir} \
            -n multiqc_report \
            {params.searchdir} \
        2> {log}
        """
