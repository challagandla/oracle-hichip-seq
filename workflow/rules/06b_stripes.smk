# Stage 06b — Architectural stripes.
#
# A stripe (flame / architectural line) is what you see when one anchor is held
# by cohesin while the extrusion complex reels in the DNA on one side: contact
# enrichment forms a continuous line running away from that anchor, rather than a
# dot. That is a different object from a loop, and a loop caller cannot find it —
# FitHiChIP tests discrete pixel pairs against a distance-decay background, so a
# stripe is at best fragmented into a row of marginal, individually
# unconvincing pixels.
#
# Added on the argument of Banecki et al. (Commun Biol 2025), who show stripes are
# a primary readout of a cohesin HiChIP experiment and that their detection is one
# of the things a good protocol buys you. This pipeline called loops, insulation
# and compartments but had no stripe caller at all.
#
# Read them with the anchor in mind. On a CTCF or cohesin anchor set a stripe is
# directly a loop-extrusion anchor. On H3K27ac the anchors are enhancers, extrusion
# is not what defines them, and stripes are correspondingly weaker — they are
# reported for completeness and should not carry an argument on their own.


rule stripenn_call:
    """Call stripes on the balanced matrix at `stripes.resolution`.

    Deliberately not run at the 5 kb loop resolution: a stripe is detected as an
    image feature (Canny edge detection over the contact map), and at 5 kb a
    HiChIP matrix is sparse enough that the edges being traced are mostly
    sampling noise. 10 kb is the coarsest resolution at which anchors remain
    resolved.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
    output:
        tsv = RESULTS / "stripes/{sample}/result_unfiltered.tsv",
        filtered = RESULTS / "stripes/{sample}/result_filtered.tsv",
    params:
        outdir = lambda wc: str(RESULTS / f"stripes/{wc.sample}"),
        res = config["stripes"]["resolution"],
        pval = config["stripes"]["max_pixel_pval"],
        canny = config["stripes"]["canny_sigma"],
        minlen = config["stripes"]["min_length"],
        maxw = config["stripes"]["max_width"],
    threads: config["threads"]["stripenn"]
    conda: "../envs/stripenn.yaml"
    log:
        RESULTS / "logs/stripes/{sample}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})

        stripenn compute \
            --cool {input.mcool}::/resolutions/{params.res} \
            --out {params.outdir}/ \
            --pvalue {params.pval} \
            --numcores {threads} \
            --canny {params.canny} \
            --minL {params.minlen} \
            --maxW {params.maxw} \
            > {log} 2>&1 || true

        # stripenn writes nothing when it finds no stripes, and a missing file is
        # indistinguishable from a crash to Snakemake. Materialise both tables with
        # a header either way, and let the QC summary report the count.
        for f in result_unfiltered.tsv result_filtered.tsv; do
            if [ ! -s {params.outdir}/$f ]; then
                printf 'chr\tpos1\tpos2\tchr2\tpos3\tpos4\tlength\twidth\ttotal\tMean\tmaxpixel\tnum\tstart\tend\tx\ty\th\tw\tmedpixel\tpvalue\n' \
                    > {params.outdir}/$f
                echo "no stripes called for {wildcards.sample}" >> {log}
            fi
        done
        """


rule stripe_summary:
    """One table across samples: stripe count, median length, anchor mark."""
    input:
        expand(RESULTS / "stripes/{sample}/result_filtered.tsv", sample=SAMPLE_IDS),
    output:
        tsv = RESULTS / "stripes/stripe_summary.tsv",
    params:
        samples = SAMPLE_IDS,
        marks = {s: SAMPLES.loc[s, "mark"] for s in SAMPLE_IDS},
        res = config["stripes"]["resolution"],
    conda: "../envs/pandas.yaml"
    log:
        RESULTS / "logs/stripes/summary.log",
    script:
        "../scripts/stripe_summary.py"
