# Stage 04 — Extract 1D reads from HiChIP and call peaks with MACS2
# HiChIP loop calling needs ChIP anchors. We split the dedup .pairs.gz into
# a BAM of 1D reads, then run MACS2 in narrow or broad mode per mark.

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
    log:
        RESULTS / "logs/pairs_to_1d/{sample}.log"
    shell:
        r"""
        pairtools split --output-sam - {input.pairsam} 2> {log} | \
            samtools sort -@ {threads} -o {output.bam} - 2>> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}
        """


def _macs2_mode(wc):
    mark = SAMPLES.loc[wc.sample, "mark"]
    return config["macs2"]["marks"].get(mark, {"mode": "narrow"})["mode"]

def _macs2_q(wc):
    mark = SAMPLES.loc[wc.sample, "mark"]
    return config["macs2"]["marks"].get(mark, {"q": 0.01})["q"]

def _macs2_broad_cutoff(wc):
    mark = SAMPLES.loc[wc.sample, "mark"]
    return config["macs2"]["marks"].get(mark, {}).get("broad_cutoff", 0.1)


rule macs2_peaks:
    input:
        bam = RESULTS / "bam_1d/{sample}.1d.bam"
    output:
        bed = RESULTS / "peaks/{sample}_peaks.bed",
        # MACS2 native outputs (narrowPeak or broadPeak) kept under peaks/raw/
        macs_out = RESULTS / "peaks/raw/{sample}_peaks_macs.done"
    params:
        gsize = config["macs2"]["genome_size"],
        # Per-rule flattened scalars — avoids dict-subscript fragility in shell
        mode = _macs2_mode,
        qval = _macs2_q,
        broad_cutoff = _macs2_broad_cutoff,
        outdir = RESULTS / "peaks/raw"
    threads: config["threads"]["macs2"]
    log:
        RESULTS / "logs/macs2/{sample}.log"
    shell:
        r"""
        mkdir -p {params.outdir}

        if [ "{params.mode}" = "broad" ]; then
            macs2 callpeak -t {input.bam} -f BAMPE -g {params.gsize} \
                --outdir {params.outdir} -n {wildcards.sample} \
                --broad --broad-cutoff {params.broad_cutoff} -q {params.qval} \
                --nomodel 2> {log}
            cut -f1-3 {params.outdir}/{wildcards.sample}_peaks.broadPeak | \
                sort -k1,1 -k2,2n > {output.bed}
        else
            macs2 callpeak -t {input.bam} -f BAMPE -g {params.gsize} \
                --outdir {params.outdir} -n {wildcards.sample} -q {params.qval} \
                --nomodel 2> {log}
            cut -f1-3 {params.outdir}/{wildcards.sample}_peaks.narrowPeak | \
                sort -k1,1 -k2,2n > {output.bed}
        fi

        touch {output.macs_out}
        """
