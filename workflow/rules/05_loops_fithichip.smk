# Stage 05 — FitHiChIP loop calling (Peak-to-ALL by default)
# FitHiChIP correctly models 1D ChIP bias via spline regression. Run per-sample
# at the configured bin size.

rule pairs_to_validpairs:
    """Convert pairtools .pairs.gz to HiC-Pro-style .allValidPairs for FitHiChIP."""
    input:
        pairs = RESULTS / "pairs/{sample}.dedup.pairs.gz"
    output:
        vpairs = RESULTS / "pairs/{sample}.allValidPairs"
    conda: "../envs/coreutils.yaml"
    log:
        RESULTS / "logs/pairs_to_validpairs/{sample}.log"
    shell:
        r"""
        # pairtools .pairs columns with read IDs retained are:
        # readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type ...
        # FitHiChIP needs the first seven HiC-Pro-like validPairs columns:
        # readID chr1 pos1 strand1 chr2 pos2 strand2
        zcat {input.pairs} | awk 'BEGIN{{OFS="\t"}} \
            !/^#/ {{ print $1, $2, $3, $6, $4, $5, $7 }}' \
            > {output.vpairs} 2> {log}
        """

rule fithichip_config:
    """
    Build a per-sample FitHiChIP config text file. FitHiChIP uses historical
    numeric codes, so those are explicit in config.yaml and written here.
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
        int_code = config["fithichip"]["int_type_code"],
        bias_code = config["fithichip"]["bias_type_code"],
        use_p2p = config["fithichip"].get("use_p2p_background", 1),
        itype = config["fithichip"].get("interaction_type", "Peak-to-ALL"),
        bgtype = config["fithichip"].get("background_type", "Coverage_Bias"),
        outdir = lambda wc: RESULTS / f"loops/{wc.sample}"
    run:
        outdir = Path(params.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        text = f"""
# Auto-generated FitHiChIP config for sample {wildcards.sample}
# interaction_type={params.itype}; background_type={params.bgtype}
ValidPairs={input.pairs}
PeakFile={input.peaks}
ChrSizeFile={input.chromsizes}
OutDir={outdir}/
IntType={params.int_code}
BINSIZE={params.bin_size}
LowDistThr={params.lower}
UppDistThr={params.upper}
QVALUE={params.fdr}
HiCProBasedir=
PREFIX={wildcards.sample}
Draw=1
TimeProf=1
OverWrite=1
UseP2PBackgrnd={params.use_p2p}
BiasType={params.bias_code}
""".lstrip()
        Path(output.cfg).write_text(text)

rule fithichip_run:
    """
    Run FitHiChIP. Produces the canonical interactions BED at the configured
    FDR threshold.
    """
    input:
        cfg = RESULTS / "loops/{sample}/fithichip.config"
    output:
        loops = RESULTS / f"loops/{{sample}}/{{sample}}.interactions_FitHiC_{FITHICHIP_Q_LABEL}.bed"
    threads: config["threads"]["fithichip"]
    conda: "../envs/fithichip.yaml"
    log:
        RESULTS / "logs/fithichip/{sample}.log"
    params:
        q_label = FITHICHIP_Q_LABEL
    shell:
        r"""
        # FitHiChIP ships two front-ends depending on installation method:
        #   bioconda package  →  `fithichip --cfg <file>`
        #   legacy shell      →  `FitHiChIP_HiCPro.sh -C <file>`
        # We prefer the bioconda entrypoint; fall back to the shell script.
        if command -v fithichip >/dev/null 2>&1; then
            fithichip --cfg {input.cfg} 2> {log}
        elif command -v FitHiChIP_HiCPro.sh >/dev/null 2>&1; then
            FitHiChIP_HiCPro.sh -C {input.cfg} 2> {log}
        else
            echo "ERROR: neither 'fithichip' nor 'FitHiChIP_HiCPro.sh' found in PATH." >&2
            echo "Install with: mamba install -c bioconda fithichip" >&2
            exit 1
        fi

        # Normalise output path — FitHiChIP versions differ in subdirectory nesting
        if [ ! -s {output.loops} ]; then
            found=$(find "$(dirname {output.loops})" -type f \
                    -name "*interactions_FitHiC_{params.q_label}.bed" | head -n 1)
            if [ -n "$found" ]; then
                cp "$found" {output.loops}
            else
                echo "ERROR: FitHiChIP produced no output BED at q={params.q_label}" >&2
                exit 1
            fi
        fi
        test -s {output.loops}
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
    conda: "../envs/mustache.yaml"
    log:
        RESULTS / "logs/mustache/{sample}.log"
    shell:
        r"""
        mustache -f {input.mcool} -r {params.res} -p {threads} -o {output.tsv} 2> {log}
        """
