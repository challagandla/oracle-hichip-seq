# Stage 03 — Build .cool / .mcool matrices and ICE-balance them
# We bin at multiple resolutions; ORACLE consumes 5 kb / 25 kb / 100 kb / 1 Mb.

BIN_SIZES_BP = [int(k) * 1000 for k in config["cooler"]["bin_sizes_kb"]]
BASE_BIN_BP = min(BIN_SIZES_BP)

rule cooler_cload:
    """
    Load pairs at the finest resolution. All coarser resolutions are derived
    via cooler zoomify (consistent + reproducible).
    """
    input:
        pairs = RESULTS / "pairs/{sample}.dedup.pairs.gz",
        index = RESULTS / "pairs/{sample}.dedup.pairs.gz.px2",
        chromsizes = GENOME["chromsizes"]
    output:
        cool = temp(RESULTS / "cool/{sample}.unfiltered.cool")
    params:
        binsize = BASE_BIN_BP,
        assembly = ASSEMBLY
    threads: config["threads"]["cooler"]
    conda: "../envs/cooler.yaml"
    log:
        RESULTS / "logs/cooler_cload/{sample}.log"
    shell:
        r"""
        cooler cload pairix \
            --nproc {threads} \
            --assembly {params.assembly} \
            {input.chromsizes}:{params.binsize} \
            {input.pairs} \
            {output.cool} 2> {log}
        """


rule cooler_blacklist_filter:
    """Remove contacts touching assembly-blacklisted bins before any analysis."""
    input:
        cool = RESULTS / "cool/{sample}.unfiltered.cool",
        blacklist = GENOME["blacklist"],
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        cool = RESULTS / "cool/{sample}.base.cool",
        json = RESULTS / "qc/blacklist/{sample}.blacklist_filter.json",
    params:
        assembly = ASSEMBLY,
    conda: "../envs/coolerpy.yaml"
    log:
        RESULTS / "logs/cooler_blacklist/{sample}.log",
    script:
        "../scripts/filter_cool_blacklist.py"


rule cooler_zoomify:
    """
    Zoomify to multi-resolution .mcool covering all configured bin sizes.
    Each resolution is ICE-balanced, its HDF5 convergence attributes are
    audited, and only proven-converged weights remain in the published mcool.
    """
    input:
        cool = RESULTS / "cool/{sample}.base.cool",
        blacklist = GENOME["blacklist"],
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        mcool = RESULTS / "cool/{sample}.mcool",
        json = RESULTS / "qc/balance/{sample}.balance.json",
        tsv = RESULTS / "qc/balance/{sample}.balance.tsv",
    params:
        resolutions_bp = sorted(BIN_SIZES_BP),
        weight_name = config["cooler"]["balance"].get("weight_name", "weight"),
        ignore_diags = config["cooler"]["balance"]["ignore_diags"],
        min_nnz = config["cooler"]["balance"]["min_nnz"],
        mad_max = config["cooler"]["balance"]["mad_max"],
        tolerance = config["cooler"]["balance"].get("tolerance", 1e-5),
        max_iterations = config["cooler"]["balance"].get("max_iterations", 200),
    threads: config["threads"]["cooler"]
    conda: "../envs/cooler.yaml"
    log:
        RESULTS / "logs/cooler_zoomify/{sample}.log"
    script:
        "../scripts/build_mcool.py"
