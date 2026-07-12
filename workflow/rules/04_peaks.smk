# Stage 04 — Extract 1D reads from HiChIP and call anchors with MACS3
# HiChIP loop calling needs ChIP anchors. We split the dedup .pairs.gz into
# a BAM of 1D reads, then call anchors in narrow or broad mode per mark.
#
# The anchors come from the HiChIP library itself, so no companion ChIP-seq
# experiment is required (Banecki et al., Commun Biol 2025, take the same
# approach in nf-HiChIP).
#
# MACS3, not MACS2: see envs/macs3.yaml. MACS2's bioconda build no longer
# loads on current glibc.
#
# Anchors are called with -f BAM, NOT -f BAMPE, and this is the substantive point
# of the stage. The two ends of a HiChIP pair are ligation partners — they sit
# megabases apart, or on different chromosomes, because that is the contact the
# assay is measuring. They are not the two ends of one sonicated fragment. Asking
# MACS to read them as a pair asks it to treat the ligation distance as a fragment
# length; MACS3 refuses outright (ZeroDivisionError building the PE track, having
# found no valid fragment). Each read end is therefore counted on its own and
# extended to a nucleosome, which is what FitHiChIP's own peak inference does
# (`--nomodel --extsize 147`).

rule pairs_to_1d_bam:
    """
    Project the DEDUPED UU pairsam to a 1D BAM (both ends, sorted + indexed).
    HiChIP peaks must come from deduped, uniquely-mapped reads only.
    """
    input:
        # Read from the deduped pairsam written by stage 02 (post-dedup, SAM intact)
        pairsam = RESULTS / "pairs/{sample}.dedup.pairsam.gz"
    output:
        bam = RESULTS / "bam_1d/{sample}.1d.bam",
        bai = RESULTS / "bam_1d/{sample}.1d.bam.bai"
    threads: config["threads"]["pairtools"]
    params:
        keep_expr = PAIRTOOLS_KEEP_EXPR
    conda: "../envs/align.yaml"
    log:
        RESULTS / "logs/pairs_to_1d/{sample}.log"
    shell:
        r"""
        pairtools select '{params.keep_expr}' \
            --output-rest /dev/null \
            --output - \
            {input.pairsam} 2> {log} | \
        pairtools split --output-sam - - 2>> {log} | \
            samtools sort -@ {threads} -o {output.bam} - 2>> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}
        """


def _macs3_mode(wc):
    mark = SAMPLES.loc[wc.sample, "mark"]
    return config["macs3"]["marks"].get(mark, {"mode": "narrow"})["mode"]

def _macs3_q(wc):
    mark = SAMPLES.loc[wc.sample, "mark"]
    return config["macs3"]["marks"].get(mark, {"q": 0.01})["q"]

def _macs3_broad_cutoff(wc):
    mark = SAMPLES.loc[wc.sample, "mark"]
    return config["macs3"]["marks"].get(mark, {}).get("broad_cutoff", 0.1)


rule macs3_peaks:
    input:
        bam = RESULTS / "bam_1d/{sample}.1d.bam"
    output:
        bed = RESULTS / "peaks/{sample}_peaks.bed",
        # MACS2 native outputs (narrowPeak or broadPeak) kept under peaks/raw/
        macs_out = RESULTS / "peaks/raw/{sample}_peaks_macs.done"
    params:
        gsize = config["macs3"]["genome_size"],
        extsize = config["macs3"].get("extsize", 147),
        # Per-rule flattened scalars — avoids dict-subscript fragility in shell
        mode = _macs3_mode,
        qval = _macs3_q,
        broad_cutoff = _macs3_broad_cutoff,
        outdir = RESULTS / "peaks/raw"
    threads: config["threads"]["macs3"]
    conda: "../envs/macs3.yaml"
    log:
        RESULTS / "logs/macs3/{sample}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}

        # Read ends as a BED of individual tags. MACS3's own BAM reader returns
        # zero tags from this BAM (it then dies dividing by the tag count), and
        # -f BED is in any case the input FitHiChIP's peak inference uses, so the
        # anchors here and the anchors FitHiChIP would infer are defined the same
        # way.
        reads={params.outdir}/{wildcards.sample}.reads.bed
        bedtools bamtobed -i {input.bam} > $reads 2> {log}

        if [ "{params.mode}" = "broad" ]; then
            macs3 callpeak -t $reads -f BED -g {params.gsize} \
                --outdir {params.outdir} -n {wildcards.sample} \
                --broad --broad-cutoff {params.broad_cutoff} -q {params.qval} \
                --nomodel --extsize {params.extsize} 2>> {log}
            cut -f1-3 {params.outdir}/{wildcards.sample}_peaks.broadPeak | \
                sort -k1,1 -k2,2n > {output.bed}
        else
            macs3 callpeak -t $reads -f BED -g {params.gsize} \
                --outdir {params.outdir} -n {wildcards.sample} -q {params.qval} \
                --nomodel --extsize {params.extsize} 2>> {log}
            cut -f1-3 {params.outdir}/{wildcards.sample}_peaks.narrowPeak | \
                sort -k1,1 -k2,2n > {output.bed}
        fi

        rm -f $reads
        test -s {output.bed}
        echo "{wildcards.sample}: $(wc -l < {output.bed}) anchors" >> {log}
        touch {output.macs_out}
        """
