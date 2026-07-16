# Stage 08 — Publication-grade visualisation
# pyGenomeTracks composite plots (arcs + heatmap + tracks), APA, virtual 4C.

rule pygenometracks_region:
    """
    Composite figure for a region of interest: matrix heatmap +
    exploratory local-insulation track + 1D peaks + FitHiChIP loop arcs +
    GENCODE gene models. Insulation from a mark-enriched HiChIP matrix is a
    local visual aid, not an unbiased TAD-boundary call.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        loops = RESULTS / f"loops/{{sample}}/{{sample}}.interactions_FitHiC_{FITHICHIP_Q_LABEL}.bed",
        peaks = RESULTS / "peaks/{sample}_peaks.bed",
        insul = RESULTS / "qc/insulation/{sample}.insulation.tsv",
        balance = RESULTS / "qc/balance/{sample}.balance.json",
        gtf = GENOME["gtf"],
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        ini = RESULTS / "viz/tracks/{sample}_{region}.ini",
        insulation_bedgraph = RESULTS / "viz/tracks/{sample}_{region}.insulation.bdg",
        loop_bedpe = RESULTS / "viz/tracks/{sample}_{region}.loops.bedpe",
        png = RESULTS / "viz/{sample}_{region}.png"
    params:
        region = lambda wc: next(r for r in config["viz"]["regions"] if r["name"] == wc.region),
        mark = sample_mark,
        res = 10000
    threads: 2
    conda: "../envs/pygenometracks.yaml"
    log:
        RESULTS / "logs/viz/{sample}_{region}.log"
    script:
        "../scripts/pygenometracks_loops.py"


rule virtual_4c:
    """
    Virtual 4C from an explicit viewpoint, providing a one-dimensional view of
    contact enrichment across the selected locus.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        balance = RESULTS / "qc/balance/{sample}.balance.json",
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        bw = RESULTS / "viz/virtual_4c/{sample}_{region}.v4c.bw",
        png = RESULTS / "viz/virtual_4c/{sample}_{region}.v4c.png"
    params:
        region = lambda wc: next(r for r in config["viz"]["regions"] if r["name"] == wc.region),
        res = 5000
    threads: 2
    conda: "../envs/coolerpy.yaml"
    log:
        RESULTS / "logs/virtual_4c/{sample}_{region}.log"
    script:
        "../scripts/virtual_4c.py"
