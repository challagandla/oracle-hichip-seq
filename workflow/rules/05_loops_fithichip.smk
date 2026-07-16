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
FITHICHIP_SHA256 = "0ab11a130ad6b070f82ce7dff330f877c4b0156f9fe90b6e582e462339380e6c"
FITHICHIP_DIR = RESULTS / f"tools/FitHiChIP-{FITHICHIP_VERSION}"


rule fithichip_install:
    """Fetch the pinned FitHiChIP release.

    Pinned to a tag, not master: FitHiChIP's config keys and output directory
    nesting have both changed between versions, and this workflow writes the
    former and parses the latter.
    """
    output:
        script = FITHICHIP_DIR / "FitHiChIP_HiCPro.sh",
        manifest = FITHICHIP_DIR / ".oracle-install.sha256",
    params:
        url = f"https://github.com/ay-lab/FitHiChIP/archive/refs/tags/{FITHICHIP_VERSION}.tar.gz",
        sha256 = FITHICHIP_SHA256,
        version = FITHICHIP_VERSION,
        cache = f".cache/downloads/FitHiChIP-{FITHICHIP_VERSION}.tar.gz",
        dest = lambda wc, output: str(Path(output.script).parent),
    conda: "../envs/fithichip.yaml"
    log:
        RESULTS / "logs/fithichip/install.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {params.dest}) $(dirname {log})
        tmp=$(mktemp -d "$(dirname {params.dest})/.FitHiChIP-{params.version}.XXXXXX")
        trap 'rm -rf "$tmp"' EXIT
        if [ -s {params.cache} ]; then
            cp {params.cache} "$tmp/src.tar.gz"
        else
            curl -L --fail --retry 3 "{params.url}" -o "$tmp/src.tar.gz" > {log} 2>&1
        fi
        echo "{params.sha256}  $tmp/src.tar.gz" | sha256sum -c - >> {log} 2>&1
        mkdir "$tmp/source"
        tar -xzf "$tmp/src.tar.gz" -C "$tmp/source" --strip-components=1 >> {log} 2>&1
        test -s "$tmp/source/FitHiChIP_HiCPro.sh"
        (cd "$tmp/source" && \
            find FitHiChIP_HiCPro.sh src Analysis Imp_Scripts -type f -print0 | \
            sort -z | xargs -0 sha256sum > .oracle-install.sha256)
        chmod +x "$tmp/source/FitHiChIP_HiCPro.sh"
        rm -rf {params.dest}
        mv "$tmp/source" {params.dest}
        test -s {output.script}
        test -s {output.manifest}
        echo "FitHiChIP installed to {params.dest}" >> {log}
        """

# FitHiChIP's ValidPairs path requires a HiC-Pro installation to build its matrix.
# This workflow instead hands over the matrix already built with cooler:
# `COOL=` takes the 5 kb single-resolution `{sample}.base.cool`, and FitHiChIP reads
# it with `cooler dump -t pixels --join`. It must be the plain .cool file, not an
# `.mcool::/resolutions/5000` URI -- FitHiChIP validates the path with `[ ! -f ... ]`,
# which a URI fails.


rule fithichip_config:
    """
    Build a per-sample FitHiChIP config text file. FitHiChIP uses historical
    numeric codes; the workflow derives them from the readable analysis choices.
    """
    input:
        cool = RESULTS / "cool/{sample}.base.cool",
        peaks = lambda wc: RESULTS / f"peaks/consensus/{SAMPLES.loc[wc.sample, 'anchor_group']}.consensus.bed",
        chromsizes = GENOME["chromsizes"]
    output:
        cfg = RESULTS / "loops/{sample}/fithichip.config"
    params:
        bin_size = config["fithichip"]["bin_size"],
        lower = config["fithichip"]["lower_distance"],
        upper = config["fithichip"]["upper_distance"],
        fdr   = config["fithichip"]["fdr_threshold"],
        int_code = FITHICHIP_INT_CODE,
        bias_code = FITHICHIP_BIAS_CODE,
        merge_int = 1 if FITHICHIP_MERGE else 0,
        use_p2p = config["fithichip"].get("use_p2p_background", 1),
        itype = config["fithichip"].get("interaction_type", "Peak-to-ALL"),
        bgtype = config["fithichip"].get("background_type", "Coverage_Bias"),
        outdir = lambda wc: RESULTS / f"loops/{wc.sample}"
    run:
        outdir = Path(params.outdir).resolve()
        outdir.mkdir(parents=True, exist_ok=True)
        # This key set follows the FitHiChIP 11.0 config contract. MergeInt ensures
        # adjacent significant bins are represented as one merged contact.
        # Resolve every path because FitHiChIP changes into its release directory
        # before running its R and Python stages.
        text = f"""
# Auto-generated FitHiChIP config for sample {wildcards.sample}
# interaction_type={params.itype}; background_type={params.bgtype}
#
# COOL is supplied because the matrix has already been built by cooler.
ValidPairs=
Interval=
Matrix=
Bed=
HIC=
COOL={Path(input.cool).resolve()}
ChrSizeFile={Path(input.chromsizes).resolve()}
PeakFile={Path(input.peaks).resolve()}
OutDir={outdir}/
CircularGenome=0
IntType={params.int_code}
BINSIZE={params.bin_size}
LowDistThr={params.lower}
UppDistThr={params.upper}
UseP2PBackgrnd={params.use_p2p}
BiasType={params.bias_code}
MergeInt={params.merge_int}
QVALUE={params.fdr}
PREFIX={wildcards.sample}
OverWrite=1
""".lstrip()
        Path(output.cfg).write_text(text)

rule fithichip_run:
    """
    Run FitHiChIP, retain its unmerged all-interaction table, and publish the
    separately filtered q-thresholded call set used for reporting and APA.
    """
    input:
        cfg = RESULTS / "loops/{sample}/fithichip.config",
        script = FITHICHIP_DIR / "FitHiChIP_HiCPro.sh",
        install_manifest = FITHICHIP_DIR / ".oracle-install.sha256",
        blacklist = GENOME["blacklist"],
        filter_script = "workflow/scripts/filter_fithichip_loops.py",
        consensus_code = "workflow/scripts/build_consensus_loops.py",
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        raw_all = RESULTS / (
            f"loops/{{sample}}/{FITHICHIP_ALL_RESULT_DIR}/"
            "{sample}.interactions_FitHiC.bed"
        ),
        loops = RESULTS / f"loops/{{sample}}/{{sample}}.interactions_FitHiC_{FITHICHIP_Q_LABEL}.bed",
        audit = RESULTS / f"loops/{{sample}}/{{sample}}.interactions_FitHiC_{FITHICHIP_Q_LABEL}.filtering.tsv",
    threads: config["threads"]["fithichip"]
    conda: "../envs/fithichip.yaml"
    log:
        RESULTS / "logs/fithichip/{sample}.log"
    params:
        q_label = FITHICHIP_Q_LABEL,
        # The exact directory and file FitHiChIP writes the configured call set into,
        # derived from the same config values written into the .config so the path we
        # read cannot drift from the run we asked for.
        result_dir = FITHICHIP_RESULT_DIR,
        all_result_dir = FITHICHIP_ALL_RESULT_DIR,
        # Formatted here, not left as a literal: Snakemake does not expand wildcards
        # inside a params string.
        result_file = lambda wc: FITHICHIP_RESULT_FILE.format(sample=wc.sample),
        min_reads = int(config["fithichip"].get("min_reads", 0)),
        q_threshold = float(config["fithichip"]["fdr_threshold"]),
        want_desc = (
            "FitHiChIP MergeNearContacts call set: nearby significant bins are "
            "grouped, but the format can retain multiple representative rows per "
            "connected neighbourhood (fithichip.merge_nearby)."
            if FITHICHIP_MERGE else
            "Raw call set: one physical loop may appear as several adjacent bin pairs."
        ),
    shell:
        r"""
        # Isolate rule packages from user-level R and Python startup configuration.
        export R_PROFILE_USER=/dev/null
        export R_ENVIRON_USER=/dev/null
        export R_LIBS_USER=""
        export R_LIBS_SITE=""

        export PYTHONNOUSERSITE=1

        # Detect a partially deleted or modified FitHiChIP installation before a
        # long R job starts. The pinned archive hash protects installation; this
        # manifest protects the extracted scripts used at execution time.
        (cd $(dirname {input.script}) && sha256sum -c .oracle-install.sha256) \
            >> {log} 2>&1

        # FitHiChIP reports some fatal errors on stdout, so capture both streams.
        bash {input.script} -C {input.cfg} >> {log} 2>&1

        outdir="$(dirname {output.loops})"

        # FitHiChIP's official differential contract is the unmerged, unthresholded
        # PREFIX.interactions_FitHiC.bed in FitHiC_BiasCorr. This is intentionally
        # distinct from the q-filtered MergeNearContacts reporting call set below.
        raw_all="$outdir/{params.all_result_dir}/{wildcards.sample}.interactions_FitHiC.bed"
        if [ ! -s "$raw_all" ]; then
            echo "ERROR: FitHiChIP produced no unmerged all-interaction table at" >&2
            echo "       $raw_all" >&2
            exit 1
        fi

        # Take the call set from the directory the configured interaction type and
        # background actually write to. A `find -name "*interactions_FitHiC_*.bed"`
        # across the whole sample tree is NOT safe: FitHiChIP also builds an
        # ALL2ALL directory while computing the background model, so a wildcard can
        # match a different call set than the one requested -- and `find` returns
        # directory order, not sorted order, so which one it picked would vary
        # between samples.
        #
        # {params.want_desc}
        want="$outdir/{params.result_dir}/{params.result_file}"
        if [ -s "$want" ]; then
            python {input.filter_script} \
                --input "$want" \
                --blacklist {input.blacklist} \
                --min-reads {params.min_reads} \
                --q-threshold {params.q_threshold} \
                --output {output.loops} \
                --audit {output.audit}
        else
            echo "ERROR: FitHiChIP produced no BED at q={params.q_label} at" >&2
            echo "       $want" >&2
            echo "       (present: $(find "$outdir" -name '*interactions_FitHiC_*.bed' 2>/dev/null | tr '\n' ' '))" >&2
            exit 1
        fi
        test -s {output.raw_all}
        test -s {output.loops}
        """


rule normalize_fithichip_all_interactions:
    """Stream and normalize every unmerged FitHiChIP interaction without q filtering."""
    input:
        raw = RESULTS / (
            f"loops/{{sample}}/{FITHICHIP_ALL_RESULT_DIR}/"
            "{sample}.interactions_FitHiC.bed"
        ),
        blacklist = GENOME["blacklist"],
        chromsizes = GENOME["chromsizes"],
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        all_interactions = RESULTS / "loops/{sample}/{sample}.interactions_FitHiC.all.tsv.gz",
        eligible = RESULTS / "loops/{sample}/{sample}.interactions_FitHiC.eligible.tsv.gz",
        audit = RESULTS / "loops/{sample}/{sample}.interactions_FitHiC.all.audit.json",
    params:
        bin_size = int(config["fithichip"]["bin_size"]),
        lower_distance = int(config["fithichip"]["lower_distance"]),
        upper_distance = int(config["fithichip"]["upper_distance"]),
        min_count = int(config["differential"]["min_count"]),
        interaction_type = config["fithichip"]["interaction_type"],
        source_relative = lambda wc: (
            f"{FITHICHIP_ALL_RESULT_DIR}/{wc.sample}.interactions_FitHiC.bed"
        ),
    threads: 1
    conda: "../envs/pandas.yaml"
    log:
        RESULTS / "logs/fithichip/{sample}.normalize_all.log",
    script:
        "../scripts/normalize_fithichip_all.py"

rule mustache_crosscheck:
    """
    Cross-check loop calling with mustache (scale-space blob detector).
    Used for sanity; not the primary call set.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        balance = RESULTS / "qc/balance/{sample}.balance.json",
        primary = RESULTS / f"loops/{{sample}}/{{sample}}.interactions_FitHiC_{FITHICHIP_Q_LABEL}.bed",
        runtime_code = "workflow/scripts/mustache_runtime.py",
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        tsv = RESULTS / "loops/{sample}/{sample}.mustache.tsv",
        status = RESULTS / "loops/{sample}/{sample}.mustache.status.json",
    params:
        res = config.get("mustache", {}).get("resolution", 10000),
        comparison_tolerance_bins = config.get("mustache", {}).get(
            "comparison_tolerance_bins", 1
        ),
    threads: 8
    conda: "../envs/mustache.yaml"
    log:
        RESULTS / "logs/mustache/{sample}.log"
    script:
        "../scripts/mustache_balance_aware.py"
