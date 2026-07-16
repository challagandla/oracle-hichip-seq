# Stage 01 — Raw read QC and adapter/quality trimming
# Inputs:  raw paired-end FASTQ
# Outputs: FastQC HTML + fastp-trimmed FASTQs

rule fastqc_raw:
    input:
        r1 = fastq_r1,
        r2 = fastq_r2
    output:
        r1_html = RESULTS / "qc/fastqc_raw/{sample}_R1_fastqc.html",
        r2_html = RESULTS / "qc/fastqc_raw/{sample}_R2_fastqc.html",
        r1_zip = RESULTS / "qc/fastqc_raw/{sample}_R1_fastqc.zip",
        r2_zip = RESULTS / "qc/fastqc_raw/{sample}_R2_fastqc.zip",
    params:
        outdir = RESULTS / "qc/fastqc_raw",
        tmpdir = lambda wc: RESULTS / f"qc/fastqc_raw/.{wc.sample}.fastqc.tmp",
    threads: config["threads"]["fastqc"]
    conda: "../envs/qc.yaml"
    log:
        RESULTS / "logs/fastqc_raw/{sample}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}/R1 {params.tmpdir}/R2
        # Separate mate directories avoid FastQC basename collisions when users
        # keep R1/R2 in different source directories under the same filename.
        fastqc -t {threads} -o {params.tmpdir}/R1 {input.r1} 2> {log}
        fastqc -t {threads} -o {params.tmpdir}/R2 {input.r2} 2>> {log}
        # rename outputs to canonical names regardless of fastq filename
        # Derive stem robustly — strip .fastq.gz, .fq.gz, .fastq, or .fq.
        R1_STEM=$(basename {input.r1}); R1_STEM=${{R1_STEM%.fastq.gz}}; R1_STEM=${{R1_STEM%.fq.gz}}; R1_STEM=${{R1_STEM%.fastq}}; R1_STEM=${{R1_STEM%.fq}}
        R2_STEM=$(basename {input.r2}); R2_STEM=${{R2_STEM%.fastq.gz}}; R2_STEM=${{R2_STEM%.fq.gz}}; R2_STEM=${{R2_STEM%.fastq}}; R2_STEM=${{R2_STEM%.fq}}

        # Fail loudly; do not swallow output-producing step failures.
        test -f "{params.tmpdir}/R1/${{R1_STEM}}_fastqc.html" || \
            {{ echo "ERROR: FastQC did not produce ${{R1_STEM}}_fastqc.html" >&2; exit 1; }}
        test -f "{params.tmpdir}/R2/${{R2_STEM}}_fastqc.html" || \
            {{ echo "ERROR: FastQC did not produce ${{R2_STEM}}_fastqc.html" >&2; exit 1; }}
        test -f "{params.tmpdir}/R1/${{R1_STEM}}_fastqc.zip" || \
            {{ echo "ERROR: FastQC did not produce ${{R1_STEM}}_fastqc.zip" >&2; exit 1; }}
        test -f "{params.tmpdir}/R2/${{R2_STEM}}_fastqc.zip" || \
            {{ echo "ERROR: FastQC did not produce ${{R2_STEM}}_fastqc.zip" >&2; exit 1; }}

        # Publish canonical sample-based names from the sample-specific temporary
        # directory. FastQC derives names from input basenames, which need not be
        # globally unique across a user cohort.
        mv "{params.tmpdir}/R1/${{R1_STEM}}_fastqc.html" {output.r1_html}
        mv "{params.tmpdir}/R2/${{R2_STEM}}_fastqc.html" {output.r2_html}
        mv "{params.tmpdir}/R1/${{R1_STEM}}_fastqc.zip" {output.r1_zip}
        mv "{params.tmpdir}/R2/${{R2_STEM}}_fastqc.zip" {output.r2_zip}
        rm -rf {params.tmpdir}
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
