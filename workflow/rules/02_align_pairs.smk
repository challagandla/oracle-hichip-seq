# Stage 02 — Alignment + pair extraction
# bwa-mem(2) -SP5M is the canonical HiC mode (skip mate rescue, soft-clip 5'
# supplementary, mark short hits as secondary). pairtools converts the BAM
# into a .pairs.gz with the canonical pair type. UU == both ends uniquely
# mapped (the trustworthy class for HiChIP analysis).

BWA_BIN = "bwa-mem2" if config["bwa"]["use_bwamem2"] else "bwa"
BWA_IDX = GENOME["bwamem2_index"] if config["bwa"]["use_bwamem2"] else GENOME["bwa_index"]
BWA_INDEX_FILES = (
    [f"{BWA_IDX}{suffix}" for suffix in (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")]
    if config["bwa"]["use_bwamem2"]
    else [f"{BWA_IDX}{suffix}" for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa")]
)
PAIRTOOLS_KEEP_EXPR = " or ".join(f'pair_type=="{pt}"' for pt in config["pairtools"]["keep_pair_types"])
PAIRTOOLS_VALID_LIGATION_EXPR = (
    "(rfrag1 >= 0) and (rfrag2 >= 0) and "
    "((chrom1 != chrom2) or (abs(rfrag1 - rfrag2) > 1))"
    if config["pairtools"].get("filter_restriction_artifacts", True)
    else "True"
)

rule bwa_align_sort_pairs:
    """
    Streamed: bwa-mem(2) -SP5M | pairtools parse → sort by chrom/pos.
    Read IDs are intentionally retained so FitHiChIP validPairs conversion can
    emit a valid HiC-Pro-like first column downstream.
    """
    input:
        r1 = RESULTS / "trimmed/{sample}_R1.trim.fastq.gz",
        r2 = RESULTS / "trimmed/{sample}_R2.trim.fastq.gz",
        chromsizes = GENOME["chromsizes"],
        digest = GENOME["digest_bed"],
        index = BWA_INDEX_FILES,
    output:
        # This large sorted intermediate is consumed only by pairtools_dedup.
        pairsam = temp(RESULTS / "pairs/{sample}.sorted.pairsam.gz")
    params:
        idx = BWA_IDX,
        bwa = BWA_BIN,
        flags = config["bwa"]["flags"],
        min_mapq = config["pairtools"]["min_mapq"],
        walks = config["pairtools"]["walks_policy"]
    threads: config["threads"]["bwa"]
    conda: "../envs/align.yaml"
    log:
        RESULTS / "logs/bwa_pairs/{sample}.log"
    shell:
        r"""
        ({params.bwa} mem {params.flags} -t {threads} {params.idx} {input.r1} {input.r2} | \
          pairtools parse \
            --chroms-path {input.chromsizes} \
            --min-mapq {params.min_mapq} \
            --walks-policy {params.walks} \
            --add-columns mapq | \
          # Annotate the MboI fragment at each end. These columns support explicit
          # dangling-end/self-circle QC while leaving the contact set unchanged.
          pairtools restrict --frags {input.digest} | \
          # SAM kept so pairtools split can emit a 1D BAM for MACS3 in stage 04.
          pairtools sort --nproc {threads} \
              --tmpdir "${{TMPDIR:-/tmp}}" \
              -o {output.pairsam}) \
          2> {log}
        """

rule pairtools_dedup:
    """
    Deduplicate at pair-level (NOT read-level). Picard MarkDuplicates is wrong
    for HiC/HiChIP; always use pairtools dedup.
    Emits restriction-filtered UU contacts, their pairix index, and duplicate
    statistics with explicit denominators for QC.
    """
    input:
        pairsam = RESULTS / "pairs/{sample}.sorted.pairsam.gz"
    output:
        pairs = RESULTS / "pairs/{sample}.dedup.pairs.gz",
        index = RESULTS / "pairs/{sample}.dedup.pairs.gz.px2",
        # Pre-filter UU contacts are retained until restriction QC completes.
        all_pairs = temp(RESULTS / "pairs/{sample}.dedup.UU.pairs.gz"),
        # This pairsam is consumed only by the 1D BAM extraction stage.
        pairsam_dedup = temp(RESULTS / "pairs/{sample}.dedup.pairsam.gz"),
        stats = RESULTS / "qc/pairtools/{sample}.dedup.stats.txt",
        # Only the unmapped count is retained in the statistics report.
        unmapped = temp(RESULTS / "pairs/{sample}.unmapped.pairs.gz")
    threads: config["threads"]["pairtools"]
    params:
        keep_expr = PAIRTOOLS_KEEP_EXPR,
        valid_ligation_expr = PAIRTOOLS_VALID_LIGATION_EXPR,
    conda: "../envs/align.yaml"
    log:
        RESULTS / "logs/pairtools_dedup/{sample}.log"
    shell:
        r"""
        # Step 1: deduplicate at pair-level (NOT read-level)
        pairtools dedup \
            --mark-dups \
            --output-stats {output.stats} \
            --output-unmapped {output.unmapped} \
            --output {output.pairsam_dedup} \
            {input.pairsam} 2> {log}

        # Step 2: keep configured high-confidence pair types and split to .pairs.gz
        # NOTE: {output.pairsam_dedup} is used as INPUT here (already written above).
        # These are sequential shell commands — no parallelism issue.
        pairtools select '{params.keep_expr}' \
            --output-rest /dev/null \
            --output - \
            {output.pairsam_dedup} 2>> {log} | \
        pairtools split --output-pairs {output.all_pairs} - 2>> {log}

        # Community-default contact maps exclude neighbouring-fragment dangling
        # ends, self-circles, same-strand mirror pairs, and unassigned fragments.
        # The unfiltered UU set remains available to restriction QC and 1D peak
        # calling still uses the deduplicated pairsam above.
        pairtools select \
            --type-cast rfrag1 int --type-cast rfrag2 int \
            --cmd-out "bgzip -c -@ {threads}" \
            --output-rest /dev/null --output {output.pairs} \
            '{params.valid_ligation_expr}' {output.all_pairs} 2>> {log}

        # Step 3: index the explicitly BGZF-compressed valid-ligation pairs. Do
        # not rely on pairtools' compressor auto-detection: some pbgzip builds
        # emit ordinary gzip, which pairix cannot index.
        test -s {output.pairs}
        pairix -f -p pairs {output.pairs} 2>> {log}
        test -s {output.index}
        """

rule pairtools_stats:
    """
    Per-sample selected-contact statistics: cis/trans and distance categories.
    Raw valid-pair yield and duplicate rate are assembled later from their exact
    fastp and pairtools-dedup denominator populations.
    """
    input:
        pairs = RESULTS / "pairs/{sample}.dedup.pairs.gz"
    output:
        stats = RESULTS / "qc/pairtools/{sample}.pairs.stats.txt"
    threads: 1
    conda: "../envs/align.yaml"
    log:
        RESULTS / "logs/pairtools_stats/{sample}.log"
    shell:
        r"""
        pairtools stats {input.pairs} -o {output.stats} 2> {log}
        """


rule restriction_fragment_qc:
    """Report pre-filter fragment classes and retained valid-ligation yield."""
    input:
        pairs = RESULTS / "pairs/{sample}.dedup.UU.pairs.gz",
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        json = RESULTS / "qc/restriction/{sample}.restriction.json",
    conda: "../envs/pandas.yaml"
    log:
        RESULTS / "logs/restriction_qc/{sample}.log",
    script:
        "../scripts/restriction_fragment_qc.py"
