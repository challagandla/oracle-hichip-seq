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
# Read them with the anchor in mind. On a CTCF or cohesin anchor set a stripe is
# directly a loop-extrusion anchor. On H3K27ac the anchors are enhancers, extrusion
# is not what defines them, and stripes are correspondingly weaker — they are
# reported for completeness and should not carry an argument on their own.


rule stripenn_call:
    """Call stripes on the configured raw-count matrix at `stripes.resolution`.

    Deliberately not run at the 5 kb loop resolution: a stripe is detected as an
    image feature (Canny edge detection over the contact map), and at 5 kb a
    HiChIP matrix is sparse enough that the edges being traced are mostly
    sampling noise. 10 kb is the coarsest resolution at which anchors remain
    resolved.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        view = RESULTS / "qc/view_main_chroms.bed",
    output:
        tsv = RESULTS / "stripes/{sample}/result_unfiltered.tsv",
        filtered = RESULTS / "stripes/{sample}/result_filtered.tsv",
    params:
        outdir = lambda wc: str(RESULTS / f"stripes/{wc.sample}"),
        res = config["stripes"]["resolution"],
        pval = config["stripes"]["max_pixel_pval"],
        canny = config["stripes"]["canny_sigma"],
        norm = config["stripes"].get("norm", "weight"),
        # bp -> bins. stripenn's --minL/--maxW are counts of bins, not base pairs;
        # config states them in bp so they stay meaningful if the resolution changes.
        minlen = max(1, config["stripes"]["min_length"] // config["stripes"]["resolution"]),
        maxw = max(1, config["stripes"]["max_width"] // config["stripes"]["resolution"]),
    threads: config["threads"]["stripenn"]
    conda: "../envs/stripenn.yaml"
    log:
        RESULTS / "logs/stripes/{sample}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})

        # Prefer the declared environment's shared libraries over host libraries.
        export LD_LIBRARY_PATH="${{CONDA_PREFIX}}/lib${{LD_LIBRARY_PATH:+:${{LD_LIBRARY_PATH}}}}"

        # --chrom: stripenn takes a quantile of the non-zero pixels of every
        # chromosome it is given --
        #
        #     mat = self.unbalLib.fetch(CHROM); np.quantile(mat[mat>0], quantile)
        #
        # -- and on an unplaced scaffold with no contacts that array is empty, which
        # is fatal (IndexError: index -1 is out of bounds for axis 0 with size 0).
        # Its own scaffold filter only drops names containing JH5/GL4/RANDOM, which
        # does not match hg38's GL000008.2 / KI270*, so all ~160 get through. Reuse
        # the same assembled-chromosome view that cooltools and HiCRep take; stripenn
        # already excludes chrM and chrY itself.
        CHROMS=$(cut -f1 {input.view} | paste -sd,)

        # Stripenn asks before clearing an existing output directory. Process
        # substitution answers without introducing a pipefail/SIGPIPE failure.
        #
        # --norm {params.norm}: RAW counts, deliberately. Feeding stripenn the
        # ICE-balanced matrices contain NaNs at unmappable bins, which suppress
        # Canny edge detection in this image-based method. See config.yaml.
        #
        # Preserve the caller's exit status; empty output and failure are distinct.
        stripenn compute \
            --cool {input.mcool}::/resolutions/{params.res} \
            --out {params.outdir}/ \
            --chrom "$CHROMS" \
            --norm {params.norm} \
            --pvalue {params.pval} \
            --numcores {threads} \
            --canny {params.canny} \
            --minL {params.minlen} \
            --maxW {params.maxw} \
            < <(printf 'Y\n') \
            > {log} 2>&1

        # Normalize an empty successful result to the documented output schema.
        for f in result_unfiltered.tsv result_filtered.tsv; do
            if [ ! -s {params.outdir}/$f ]; then
                printf 'chr\tpos1\tpos2\tchr2\tpos3\tpos4\tlength\twidth\tMean\tmaxpixel\tpvalue\tStripiness\n' \
                    > {params.outdir}/$f
                echo "stripenn exited 0 but called no stripes for {wildcards.sample}" >> {log}
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
