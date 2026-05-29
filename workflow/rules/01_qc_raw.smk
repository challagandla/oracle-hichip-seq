# Stage 01 — Raw read QC and adapter/quality trimming
# Inputs:  raw paired-end FASTQ
# Outputs: FastQC HTML + fastp-trimmed FASTQs

rule fastqc_raw:
    input:
        r1 = fastq_r1,
        r2 = fastq_r2
    output:
        r1_html = RESULTS / "qc/fastqc_raw/{sample}_R1_fastqc.html",
        r2_html = RESULTS / "qc/fastqc_raw/{sample}_R2_fastqc.html"
    params:
        outdir = RESULTS / "qc/fastqc_raw"
    threads: config["threads"]["fastqc"]
    log:
        RESULTS / "logs/fastqc_raw/{sample}.log"
    shell:
        r"""
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2} 2> {log}
        # rename outputs to canonical names regardless of fastq filename
        mv {params.outdir}/$(basename {input.r1} .fastq.gz)_fastqc.html {output.r1_html} 2>>{log} || true
        mv {params.outdir}/$(basename {input.r2} .fastq.gz)_fastqc.html {output.r2_html} 2>>{log} || true
        """


rule fastp_trim:
    """
    Adapter + quality trim with fastp. HiChIP libraries are paired-end and
    often have adapter readthrough on short fragments — never skip.
    """
    input:
        r1 = fastq_r1,
        r2 = fastq_r2
    output:
        r1 = RESULTS / "trimmed/{sample}_R1.trim.fastq.gz",
        r2 = RESULTS / "trimmed/{sample}_R2.trim.fastq.gz",
        json = RESULTS / "qc/fastp/{sample}.fastp.json",
        html = RESULTS / "qc/fastp/{sample}.fastp.html"
    threads: config["threads"]["fastp"]
    params:
        q   = config["fastp"]["quality_cutoff"],
        l   = config["fastp"]["length_required"],
        a1  = config["fastp"]["adapter_sequence"],
        a2  = config["fastp"]["adapter_sequence_r2"]
    log:
        RESULTS / "logs/fastp/{sample}.log"
    shell:
        r"""
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              --adapter_sequence {params.a1} --adapter_sequence_r2 {params.a2} \
              --qualified_quality_phred {params.q} --length_required {params.l} \
              --detect_adapter_for_pe \
              --json {output.json} --html {output.html} \
              --thread {threads} \
              2> {log}
        """
