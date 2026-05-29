# Stage 05 — FitHiChIP loop calling (Peak-to-ALL)
# FitHiChIP is the gold standard for HiChIP loops; it correctly models the
# 1D ChIP bias via a spline regression. Run per-sample at 5 kb bin size.
#
# Note on input format:
#   FitHiChIP historically consumes the HiC-Pro `.allValidPairs` TSV. We
#   convert from pairtools .pairs.gz with rule `pairs_to_validpairs` below.
#   The conversion drops header lines and reorders columns to the HiC-Pro
#   schema: read_id, chr1, pos1, strand1, chr2, pos2, strand2, frag_size,
#           valid_type, frag1, frag2

rule pairs_to_validpairs:
    """Convert pairtools .pairs.gz to HiC-Pro .allValidPairs format for FitHiChIP."""
    input:
        pairs = RESULTS / "pairs/{sample}.dedup.pairs.gz"
    output:
        vpairs = RESULTS / "pairs/{sample}.allValidPairs"
    log:
        RESULTS / "logs/pairs_to_validpairs/{sample}.log"
    shell:
        r"""
        # Strip pairtools header (starts with #), then emit HiC-Pro 7-col validpairs:
        #   id  chr1  pos1  strand1  chr2  pos2  strand2  (FitHiChIP only needs first 7)
        zcat {input.pairs} | awk 'BEGIN{{OFS="\t"}} \
            !/^#/ {{ print $1, $2, $3, $6, $4, $5, $7 }}' \
            > {output.vpairs} 2> {log}
        """


rule fithichip_config:
    """
    Build a per-sample FitHiChIP config text file. FitHiChIP requires a
    self-contained INI-like config; we templatise it from snakemake.
    """
    input:
        pairs = RESULTS / "pairs/{sample}.allValidPairs",
        peaks = RESULTS / "peaks/{sample}_peaks.bed",
        chromsizes = GENOME["chromsizes"]
    output:
        cfg = RESULTS / "loops/{sample}/fithichip.config"
    params:
        bin_size = config["fithichip"]["bin_size"],
        lower = config["fithichip"]["lower_distance"],
        upper = config["fithichip"]["upper_distance"],
        fdr   = config["fithichip"]["fdr_threshold"],
        itype = config["fithichip"]["interaction_type"],
        bgtype = config["fithichip"]["background_type"],
        outdir = lambda wc: RESULTS / f"loops/{wc.sample}"
    run:
        outdir = Path(params.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        text = f"""
# Auto-generated FitHiChIP config for sample {wildcards.sample}
ValidPairs={input.pairs}
PeakFile={input.peaks}
ChrSizeFile={input.chromsizes}
OutDir={outdir}/
IntType=3                       # Peak-to-ALL
BINSIZE={params.bin_size}
LowDistThr={params.lower}
UppDistThr={params.upper}
QVALUE={params.fdr}
HiCProBasedir=
PREFIX={wildcards.sample}
Draw=1
TimeProf=1
OverWrite=1
UseP2PBackgrnd=1
BiasType=1                       # coverage bias correction (recommended)
""".lstrip()
        Path(output.cfg).write_text(text)


rule fithichip_run:
    """
    Run FitHiChIP. Produces the canonical interactions BED at the
    configured FDR threshold.
    """
    input:
        cfg = RESULTS / "loops/{sample}/fithichip.config"
    output:
        loops = RESULTS / "loops/{sample}/{sample}.interactions_FitHiC_Q0.01.bed"
    threads: config["threads"]["fithichip"]
    log:
        RESULTS / "logs/fithichip/{sample}.log"
    shell:
        r"""
        # FitHiChIP_HiCPro.sh is the runnable script shipped by the bioconda package
        FitHiChIP_HiCPro.sh -C {input.cfg} 2> {log}
        """


rule mustache_crosscheck:
    """
    Cross-check loop calling with mustache (scale-space blob detector).
    Used for sanity; not the primary call set.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool"
    output:
        tsv = RESULTS / "loops/{sample}/{sample}.mustache.tsv"
    params:
        res = config["fithichip"]["bin_size"]
    threads: 8
    log:
        RESULTS / "logs/mustache/{sample}.log"
    shell:
        r"""
        mustache -f {input.mcool} -r {params.res} -p {threads} -o {output.tsv} 2> {log}
        """
