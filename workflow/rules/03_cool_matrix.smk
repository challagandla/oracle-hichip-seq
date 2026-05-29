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
        chromsizes = GENOME["chromsizes"]
    output:
        cool = RESULTS / "cool/{sample}.base.cool"
    params:
        binsize = BASE_BIN_BP
    threads: config["threads"]["cooler"]
    log:
        RESULTS / "logs/cooler_cload/{sample}.log"
    shell:
        r"""
        cooler cload pairix \
            --nproc {threads} \
            --assembly {ASSEMBLY} \
            {input.chromsizes}:{params.binsize} \
            {input.pairs} \
            {output.cool} 2> {log}
        """


rule cooler_zoomify:
    """
    Zoomify to multi-resolution .mcool covering all configured bin sizes.
    Each resolution is ICE-balanced.
    """
    input:
        cool = RESULTS / "cool/{sample}.base.cool"
    output:
        mcool = RESULTS / "cool/{sample}.mcool"
    params:
        resolutions = ",".join(str(b) for b in sorted(BIN_SIZES_BP)),
        ignore_diags = config["cooler"]["balance"]["ignore_diags"],
        min_nnz = config["cooler"]["balance"]["min_nnz"],
        mad_max = config["cooler"]["balance"]["mad_max"]
    threads: config["threads"]["cooler"]
    log:
        RESULTS / "logs/cooler_zoomify/{sample}.log"
    shell:
        r"""
        cooler zoomify \
            --nproc {threads} \
            --balance \
            --balance-args "--mad-max {params.mad_max} --min-nnz {params.min_nnz} --ignore-diags {params.ignore_diags}" \
            --resolutions {params.resolutions} \
            -o {output.mcool} \
            {input.cool} 2> {log}
        """
