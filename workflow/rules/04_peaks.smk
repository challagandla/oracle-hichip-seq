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


def _primary_chrom_regex(wc):
    """Exact assembled autosomes plus chrX for the configured organism."""
    if GENOME["macs3_genome_size"] == "mm":
        return r"^chr([1-9]|1[0-9]|X)$"
    return r"^chr([1-9]|1[0-9]|2[0-2]|X)$"


rule macs3_peaks:
    input:
        bam = RESULTS / "bam_1d/{sample}.1d.bam",
        blacklist = GENOME["blacklist"],
        # The inline decompressor imports only utils.py.  Keeping this dependency
        # exact avoids invalidating expensive peak calls when unrelated shared
        # visualization or balance helpers change.
        utils = "workflow/scripts/utils.py",
    output:
        bed = RESULTS / "peaks/{sample}_peaks.bed",
        # MACS3 native outputs (narrowPeak or broadPeak) kept under peaks/raw/
        macs_out = RESULTS / "peaks/raw/{sample}_peaks_macs.done"
    params:
        gsize = GENOME["macs3_genome_size"],
        extsize = config["macs3"].get("extsize", 147),
        # Per-rule flattened scalars — avoids dict-subscript fragility in shell
        mode = _macs3_mode,
        qval = _macs3_q,
        broad_cutoff = _macs3_broad_cutoff,
        primary_regex = _primary_chrom_regex,
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
        bedtools bamtobed -i {input.bam} 2> {log} | \
            awk -v primary='{params.primary_regex}' '$1 ~ primary' > $reads
        test -s $reads

        if [ "{params.mode}" = "broad" ]; then
            macs3 callpeak -t $reads -f BED -g {params.gsize} \
                --outdir {params.outdir} -n {wildcards.sample} \
                --broad --broad-cutoff {params.broad_cutoff} -q {params.qval} \
                --nomodel --extsize {params.extsize} --keep-dup all 2>> {log}
            cut -f1-3 {params.outdir}/{wildcards.sample}_peaks.broadPeak > \
                {output.bed}.all_contigs
        else
            macs3 callpeak -t $reads -f BED -g {params.gsize} \
                --outdir {params.outdir} -n {wildcards.sample} -q {params.qval} \
                --nomodel --extsize {params.extsize} --keep-dup all 2>> {log}
            cut -f1-3 {params.outdir}/{wildcards.sample}_peaks.narrowPeak > \
                {output.bed}.all_contigs
        fi

        # Primary assembled chromosomes only. Alternate/unplaced contigs have
        # unstable mappability and cannot support a shared biological anchor
        # universe across libraries.
        awk -v primary='{params.primary_regex}' \
            'BEGIN{{OFS="\t"}} $1 ~ primary {{print $1,$2,$3}}' \
            {output.bed}.all_contigs | sort -u -k1,1V -k2,2n \
            > {output.bed}.unfiltered

        # Blacklisted anchors otherwise seed artifactual loops and propagate into
        # the differential candidate set. Materialise by gzip magic bytes rather
        # than a filename suffix: downloaded references are frequently renamed.
        python -c 'import shutil,sys; sys.path.insert(0,"workflow/scripts"); from utils import open_text_auto; src=open_text_auto(sys.argv[1]); dst=open(sys.argv[2],"w"); shutil.copyfileobj(src,dst); src.close(); dst.close()' \
            {input.blacklist} {output.bed}.blacklist.tmp
        bedtools intersect -nonamecheck -v -a {output.bed}.unfiltered \
            -b {output.bed}.blacklist.tmp > {output.bed}

        rm -f $reads {output.bed}.all_contigs {output.bed}.unfiltered \
            {output.bed}.blacklist.tmp
        test -s {output.bed}
        echo "{wildcards.sample}: $(wc -l < {output.bed}) anchors" >> {log}
        touch {output.macs_out}
        """


def _anchor_group_peak_inputs(wc):
    samples = SAMPLES[
        SAMPLES["anchor_group"] == wc.anchor_group
    ]["sample_id"].tolist()
    return expand(RESULTS / "peaks/{sample}_peaks.bed", sample=samples)


def _anchor_group_peak_support(wc):
    return 2 if len(_anchor_group_peak_inputs(wc)) >= 2 else 1


rule consensus_peaks:
    """A common anchor universe within one configured assay stratum."""
    input:
        peaks = _anchor_group_peak_inputs,
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        bed = RESULTS / "peaks/consensus/{anchor_group}.consensus.bed",
        audit = RESULTS / "peaks/consensus/{anchor_group}.consensus_support.tsv",
    params:
        min_support = _anchor_group_peak_support,
    conda: "../envs/pandas.yaml"
    log:
        RESULTS / "logs/consensus_peaks/{anchor_group}.log",
    script:
        "../scripts/consensus_peaks.py"


rule anchor_frip_qc:
    """Read-end FRiP in sample peaks and assay-stratum consensus anchors."""
    input:
        bam = RESULTS / "bam_1d/{sample}.1d.bam",
        sample_peaks = RESULTS / "peaks/{sample}_peaks.bed",
        consensus_peaks = lambda wc: RESULTS / f"peaks/consensus/{SAMPLES.loc[wc.sample, 'anchor_group']}.consensus.bed",
    output:
        tsv = RESULTS / "qc/anchors/{sample}.anchor_qc.tsv",
    params:
        primary_regex = _primary_chrom_regex,
    conda: "../envs/macs3.yaml"
    log:
        RESULTS / "logs/anchor_qc/{sample}.log",
    shell:
        r"""
        set -euo pipefail
        primary_bam={output.tsv}.primary.tmp.bam
        trap 'rm -f "$primary_bam"' EXIT
        regions=$(samtools idxstats {input.bam} | \
            awk -v primary='{params.primary_regex}' '$1 ~ primary {{printf "%s ",$1}}')
        test -n "$regions"
        samtools view -b -F 4 {input.bam} $regions > $primary_bam
        total=$(samtools view -c $primary_bam)
        in_sample=$(bedtools intersect -nonamecheck -u -abam $primary_bam \
            -b {input.sample_peaks} | samtools view -c)
        in_consensus=$(bedtools intersect -nonamecheck -u -abam $primary_bam \
            -b {input.consensus_peaks} | samtools view -c)
        n_sample=$(wc -l < {input.sample_peaks})
        n_consensus=$(wc -l < {input.consensus_peaks})
        awk -v t="$total" -v ps="$in_sample" -v pc="$in_consensus" \
            -v ns="$n_sample" -v nc="$n_consensus" \
            'BEGIN{{OFS="\t"; print "analysis_population","total_read_ends","read_ends_in_sample_peaks","sample_peak_frip","n_sample_peaks","read_ends_in_consensus_peaks","consensus_peak_frip","n_consensus_peaks"; print "primary_autosomes_chrX",t,ps,(t ? ps/t : 0),ns,pc,(t ? pc/t : 0),nc}}' \
            > {output.tsv}
        rm -f $primary_bam
        trap - EXIT
        echo "{wildcards.sample}: sample peaks $in_sample/$total ($n_sample peaks); consensus peaks $in_consensus/$total ($n_consensus peaks)" > {log}
        """
