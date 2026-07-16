# Stage 10 — Aggregate report
# MultiQC pulls FastQC, fastp, pairtools, MACS3, FitHiChIP, HiCRep,
# APA, and ORACLE QC into a single HTML.

rule apa_multiqc_content:
    input:
        json = RESULTS / "qc/apa/{sample}.apa.json"
    output:
        json = RESULTS / "qc/apa/{sample}_apa_mqc.json"
    params:
        kind = "apa"
    conda: "../envs/pandas.yaml"
    log:
        RESULTS / "logs/multiqc_content/{sample}.apa.log"
    script:
        "../scripts/multiqc_content.py"


rule loop_qc_multiqc_content:
    input:
        json = RESULTS / "qc/loop_qc/{sample}.json"
    output:
        json = RESULTS / "qc/loop_qc/{sample}_loop_qc_mqc.json"
    params:
        kind = "loop_qc"
    conda: "../envs/pandas.yaml"
    log:
        RESULTS / "logs/multiqc_content/{sample}.loop_qc.log"
    script:
        "../scripts/multiqc_content.py"


rule balance_multiqc_content:
    input:
        json = RESULTS / "qc/balance/{sample}.balance.json"
    output:
        json = RESULTS / "qc/balance/{sample}_balance_mqc.json"
    params:
        kind = "balance"
    conda: "../envs/pandas.yaml"
    log:
        RESULTS / "logs/multiqc_content/{sample}.balance.log"
    script:
        "../scripts/multiqc_content.py"


rule differential_multiqc_content:
    input:
        json = RESULTS / "diff/{comparison}/design.json"
    output:
        json = RESULTS / "diff/{comparison}/design_mqc.json"
    params:
        kind = "differential"
    conda: "../envs/pandas.yaml"
    log:
        RESULTS / "logs/multiqc_content/{comparison}.differential.log"
    script:
        "../scripts/multiqc_content.py"


rule multiqc_report:
    input:
        mqc_cfg = "config/multiqc_config.yaml",
        fastqc = expand(
            RESULTS / "qc/fastqc_raw/{sample}_{mate}_fastqc.zip",
            sample=SAMPLE_IDS,
            mate=["R1", "R2"],
        ),
        fastp  = expand(RESULTS / "qc/fastp/{sample}.fastp.json", sample=SAMPLE_IDS),
        pairs  = expand(RESULTS / "qc/pairtools/{sample}.pairs.stats.txt", sample=SAMPLE_IDS),
        dedup  = expand(RESULTS / "qc/pairtools/{sample}.dedup.stats.txt", sample=SAMPLE_IDS),
        apa_mqc = expand(RESULTS / "qc/apa/{sample}_apa_mqc.json", sample=SAMPLE_IDS),
        loop_mqc = expand(RESULTS / "qc/loop_qc/{sample}_loop_qc_mqc.json", sample=SAMPLE_IDS),
        balance_mqc = expand(RESULTS / "qc/balance/{sample}_balance_mqc.json", sample=SAMPLE_IDS),
        differential_mqc = expand(
            RESULTS / "diff/{comparison}/design_mqc.json",
            comparison=DIFF_COMPARISON_NAMES,
        ),
    output:
        html = RESULTS / "multiqc/multiqc_report.html"
    params:
        outdir = RESULTS / "multiqc"
    conda: "../envs/multiqc.yaml"
    log:
        RESULTS / "logs/multiqc/multiqc.log"
    shell:
        r"""
        multiqc -f --strict --no-data-dir \
            --config {input.mqc_cfg} \
            -o {params.outdir} \
            -n multiqc_report \
            {input.fastqc} {input.fastp} {input.pairs} {input.dedup} \
            {input.apa_mqc} {input.loop_mqc} {input.balance_mqc} \
            {input.differential_mqc} \
        2> {log}
        """
