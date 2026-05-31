# Stage 02 — Alignment + pair extraction
# bwa-mem(2) -SP5M is the canonical HiC mode (skip mate rescue, soft-clip 5'
# supplementary, mark short hits as secondary). pairtools converts the BAM
# into a .pairs.gz with the canonical pair type. UU == both ends uniquely
# mapped (the trustworthy class for HiChIP analysis).

BWA_BIN = "bwa-mem2" if config["bwa"]["use_bwamem2"] else "bwa"
BWA_IDX = GENOME["bwamem2_index"] if config["bwa"]["use_bwamem2"] else GENOME["bwa_index"]

rule bwa_align_sort_pairs:
    """
    Streamed: bwa-mem(2) -SP5M | pairtools parse → sort by chrom/pos.
    Read IDs are intentionally retained so FitHiChIP validPairs conversion can
    emit a valid HiC-Pro-like first column downstream.
    """
    input:
        r1 = RESULTS / "trimmed/{sample}_R1.trim.fastq.gz",
        r2 = RESULTS / "trimmed/{sample}_R2.trim.fastq.gz",
        chromsizes = GENOME["chromsizes"]
    output:
        pairsam = RESULTS / "pairs/{sample}.sorted.pairsam.gz"
    params:
        idx = BWA_IDX,
        bwa = BWA_BIN,
        flags = config["bwa"]["flags"],
        min_mapq = config["pairtools"]["min_mapq"],
        walks = config["pairtools"]["walks_policy"]
    threads: config["threads"]["bwa"]
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
          # SAM kept so pairtools split can emit a 1D BAM for MACS2 in stage 04.
          pairtools sort --nproc {threads} \
              --tmpdir "${{TMPDIR:-/tmp}}" \
              -o {output.pairsam}) \
          2> {log}
        """

rule pairtools_dedup:
    """
    Deduplicate at pair-level (NOT read-level). Picard MarkDuplicates is wrong
    for HiC/HiChIP; always use pairtools dedup.
    Emits clean dedup pairs + duplicate stats for QC.
    """
    input:
        pairsam = RESULTS / "pairs/{sample}.sorted.pairsam.gz"
    output:
        pairs = RESULTS / "pairs/{sample}.dedup.pairs.gz",
        pairsam_dedup = RESULTS / "pairs/{sample}.dedup.pairsam.gz",
        stats = RESULTS / "qc/pairtools/{sample}.dedup.stats.txt",
        unmapped = RESULTS / "pairs/{sample}.unmapped.pairs.gz"
    threads: config["threads"]["pairtools"]
    params:
        keep_types = ",".join(config["pairtools"]["keep_pair_types"])
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

        # Step 2: keep only UU pairs and split to .pairs.gz
        # NOTE: {output.pairsam_dedup} is used as INPUT here (already written above).
        # These are sequential shell commands — no parallelism issue.
        pairtools select '(pair_type=="UU")' \
            --output-rest /dev/null \
            --output - \
            {output.pairsam_dedup} 2>> {log} | \
        pairtools split --output-pairs {output.pairs} - 2>> {log}

        # Step 3: index the .pairs.gz with pairix (required by cooler cload pairix)
        pairix -f -p pairs {output.pairs} 2>> {log}
        """

rule pairtools_stats:
    """
    Per-sample pair statistics: cis/trans, distance-decay categories,
    duplicate fraction, valid-pair yield.
    """
    input:
        pairs = RESULTS / "pairs/{sample}.dedup.pairs.gz"
    output:
        stats = RESULTS / "qc/pairtools/{sample}.pairs.stats.txt"
    threads: 1
    log:
        RESULTS / "logs/pairtools_stats/{sample}.log"
    shell:
        r"""
        pairtools stats {input.pairs} -o {output.stats} 2> {log}
        """
