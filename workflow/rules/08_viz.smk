# Stage 08 — Publication-grade visualisation
# pyGenomeTracks composite plots (arcs + heatmap + tracks), APA, virtual 4C.

rule pygenometracks_region:
    """
    Composite figure for a region of interest: matrix heatmap +
    insulation track + 1D peaks + arc plot of FitHiChIP loops.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        loops = RESULTS / f"loops/{{sample}}/{{sample}}.interactions_FitHiC_{FITHICHIP_Q_LABEL}.bed",
        peaks = RESULTS / "peaks/{sample}_peaks.bed",
        insul = RESULTS / "qc/insulation/{sample}.insulation.tsv"
    output:
        ini = RESULTS / "viz/tracks/{sample}_{region}.ini",
        png = RESULTS / "viz/{sample}_{region}.png"
    params:
        region = lambda wc: next(r for r in config["viz"]["regions"] if r["name"] == wc.region),
        mark = sample_mark,
        res = 10000
    threads: 2
    log:
        RESULTS / "logs/viz/{sample}_{region}.log"
    script:
        "../scripts/pygenometracks_loops.py"


rule virtual_4c:
    """
    Virtual 4C from a viewpoint (anchor of interest) — used to show that
    a loop anchor 'looks at' specific downstream targets.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool"
    output:
        bw = RESULTS / "viz/virtual_4c/{sample}_{region}.v4c.bw",
        png = RESULTS / "viz/virtual_4c/{sample}_{region}.v4c.png"
    params:
        region = lambda wc: next(r for r in config["viz"]["regions"] if r["name"] == wc.region),
        res = 5000
    threads: 2
    log:
        RESULTS / "logs/virtual_4c/{sample}_{region}.log"
    script:
        "../scripts/virtual_4c.py"
