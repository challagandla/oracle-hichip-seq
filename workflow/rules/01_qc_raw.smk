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
    conda: "../envs/qc.yaml"
    log:
        RESULTS / "logs/fastqc_raw/{sample}.log"
    shell:
        r"""
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2} 2> {log}
        # rename outputs to canonical names regardless of fastq filename
        # Derive stem robustly — strip .fastq.gz, .fq.gz, or .fastq
        R1_STEM=$(basename {input.r1}); R1_STEM=${{R1_STEM%.fastq.gz}}; R1_STEM=${{R1_STEM%.fq.gz}}; R1_STEM=${{R1_STEM%.fastq}}
        R2_STEM=$(basename {input.r2}); R2_STEM=${{R2_STEM%.fastq.gz}}; R2_STEM=${{R2_STEM%.fq.gz}}; R2_STEM=${{R2_STEM%.fastq}}

        # Fail loudly; do not swallow output-producing step failures.
        test -f "{params.outdir}/${{R1_STEM}}_fastqc.html" || \
            {{ echo "ERROR: FastQC did not produce ${{R1_STEM}}_fastqc.html" >&2; exit 1; }}
        test -f "{params.outdir}/${{R2_STEM}}_fastqc.html" || \
            {{ echo "ERROR: FastQC did not produce ${{R2_STEM}}_fastqc.html" >&2; exit 1; }}

        # Only rename when the names actually differ. The FASTQs here are already
        # named {{sample}}_R1.fastq.gz, so the "canonical" target IS the FastQC
        # output, and `mv x x` fails with "are the same file" -- taking the rule
        # down over a rename it did not need to do.
        [ "{params.outdir}/${{R1_STEM}}_fastqc.html" -ef "{output.r1_html}" ] || \
            mv "{params.outdir}/${{R1_STEM}}_fastqc.html" {output.r1_html} 2>>{log}
        [ "{params.outdir}/${{R2_STEM}}_fastqc.html" -ef "{output.r2_html}" ] || \
            mv "{params.outdir}/${{R2_STEM}}_fastqc.html" {output.r2_html} 2>>{log}
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
        r1 = temp(RESULTS / "trimmed/{sample}_R1.trim.fastq.gz"),
        r2 = temp(RESULTS / "trimmed/{sample}_R2.trim.fastq.gz"),
        json = RESULTS / "qc/fastp/{sample}.fastp.json",
        html = RESULTS / "qc/fastp/{sample}.fastp.html"
    threads: config["threads"]["fastp"]
    params:
        q   = config["fastp"]["quality_cutoff"],
        l   = config["fastp"]["length_required"],
        a1  = config["fastp"]["adapter_sequence"],
        a2  = config["fastp"]["adapter_sequence_r2"]
    conda: "../envs/qc.yaml"
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
