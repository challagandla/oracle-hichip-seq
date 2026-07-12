# Stage 05 — FitHiChIP loop calling (Peak-to-ALL by default)
# FitHiChIP models the 1D ChIP bias via spline regression, which is what makes a
# HiChIP loop call different from a Hi-C one: contact frequency between two
# anchors is inflated simply because both are enriched in the ChIP, and without
# correcting for that the top "loops" are just the strongest peaks.
#
# FitHiChIP is not distributed through conda. It is a GitHub release of shell + R
# + python scripts, so the pinned release is fetched here and its dependencies
# come from envs/fithichip.yaml.


FITHICHIP_VERSION = "11.0"
FITHICHIP_DIR = RESULTS / f"tools/FitHiChIP-{FITHICHIP_VERSION}"


rule fithichip_install:
    """Fetch the pinned FitHiChIP release.

    Pinned to a tag, not master: FitHiChIP's config keys and output directory
    nesting have both changed between versions, and this workflow writes the
    former and parses the latter.
    """
    output:
        script = FITHICHIP_DIR / "FitHiChIP_HiCPro.sh",
    params:
        url = f"https://github.com/ay-lab/FitHiChIP/archive/refs/tags/{FITHICHIP_VERSION}.tar.gz",
        dest = lambda wc, output: str(Path(output.script).parent),
    conda: "../envs/fithichip.yaml"
    log:
        RESULTS / "logs/fithichip/install.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.dest} $(dirname {log})
        curl -L --fail --retry 3 "{params.url}" -o {params.dest}/src.tar.gz > {log} 2>&1
        tar -xzf {params.dest}/src.tar.gz -C {params.dest} --strip-components=1 >> {log} 2>&1
        rm -f {params.dest}/src.tar.gz
        chmod +x {output.script}
        test -s {output.script}
        echo "FitHiChIP installed to {params.dest}" >> {log}
        """

# The `ValidPairs=` input path is a dead end and is deliberately not used.
#
# FitHiChIP 11.0 only accepts validPairs if it can build the contact matrix and bin
# interval files from them, and it does that by shelling out to HiC-Pro:
#
#   if [[ -z $InpCoolFile && -z $InpHiCFile && -z $InpInitialInteractionBedFile ]]; then
#       if [[ -z $InpBinIntervalFile || -z $InpMatrixFile ]]; then
#           HiCProExec=`which HiC-Pro`
#           if [[ -z $HiCProExec ]]; then
#               echo 'ERROR ===>>>> HiC-pro is not installed ... FitHiChIP quits !!!'
#
# HiC-Pro is not installed, is not on bioconda, and is the legacy stack this
# pipeline exists to avoid. Feeding FitHiChIP validPairs therefore made it quit
# before calling a single loop, for every sample.
#
# We already build the matrix with cooler, so it is handed over directly:
# `COOL=` takes the 5 kb single-resolution `{sample}.base.cool`, and FitHiChIP reads
# it with `cooler dump -t pixels --join`. It must be the plain .cool file, not an
# `.mcool::/resolutions/5000` URI -- FitHiChIP validates the path with `[ ! -f ... ]`,
# which a URI fails.


rule fithichip_config:
    """
    Build a per-sample FitHiChIP config text file. FitHiChIP uses historical
    numeric codes, so those are explicit in config.yaml and written here.
    """
    input:
        cool = RESULTS / "cool/{sample}.base.cool",
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
        # Key set follows FitHiChIP 11.0's own configfile. The previous template
        # carried Draw / TimeProf / HiCProBasedir, which 11.0 does not read, and
        # omitted MergeInt, which it does: with MergeInt unset, adjacent
        # significant bins are reported as separate loops instead of being merged
        # into one contact, which inflates the loop count and double-counts the
        # same interaction in the differential test downstream.
        text = f"""
# Auto-generated FitHiChIP config for sample {wildcards.sample}
# interaction_type={params.itype}; background_type={params.bgtype}
#
# COOL, not ValidPairs: the validPairs path makes FitHiChIP shell out to HiC-Pro to
# build the matrix, and quit when it is absent. cooler already built it.
ValidPairs=
Interval=
Matrix=
Bed=
HIC=
COOL={Path(input.cool).resolve()}
ChrSizeFile={input.chromsizes}
PeakFile={input.peaks}
OutDir={outdir}/
CircularGenome=0
IntType={params.int_code}
BINSIZE={params.bin_size}
LowDistThr={params.lower}
UppDistThr={params.upper}
UseP2PBackgrnd={params.use_p2p}
BiasType={params.bias_code}
MergeInt=1
QVALUE={params.fdr}
PREFIX={wildcards.sample}
OverWrite=1
""".lstrip()
        Path(output.cfg).write_text(text)

rule fithichip_run:
    """
    Run FitHiChIP. Produces the canonical interactions BED at the configured
    FDR threshold.
    """
    input:
        cfg = RESULTS / "loops/{sample}/fithichip.config",
        script = FITHICHIP_DIR / "FitHiChIP_HiCPro.sh",
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
        # Isolate the conda R from the host installation. ~/.Rprofile here runs
        # .libPaths("~/Rlibs"), which PREPENDS a library built against a different
        # R to the search path, and FitHiChIP's R steps (edgeR, ggplot2,
        # data.table, dplyr) then fail to load packages that are in fact installed
        # in this environment. The profile is sourced on every R startup, so
        # clearing R_LIBS_USER alone does not help — the profile runs after the
        # environment is read.
        export R_PROFILE_USER=/dev/null
        export R_ENVIRON_USER=/dev/null
        export R_LIBS_USER=""
        export R_LIBS_SITE=""

        # FitHiChIP is not a conda package, so there is no binary to look for on
        # PATH — it is the release fetched by fithichip_install. The previous
        # version of this rule probed for a `fithichip` executable and, failing
        # that, told the user to `mamba install -c bioconda fithichip`, which does
        # not exist on any channel.
        bash {input.script} -C {input.cfg} 2> {log}

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
