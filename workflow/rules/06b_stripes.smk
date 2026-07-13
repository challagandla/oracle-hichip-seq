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

        # Put the environment's own libraries ahead of the host's. stripenn pulls
        # the opencv-python wheel, whose cv2 .so needs a newer libstdc++ than the
        # system provides and carries no RPATH back into the env, so the loader
        # falls through to /lib/x86_64-linux-gnu and dies on a missing CXXABI
        # version. Same failure mode as a host ~/.Rprofile or ~/.local numpy
        # shadowing a conda environment: the env is correct, the search order is not.
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

        # stdin via process substitution, NOT `yes Y | stripenn`. stripenn asks on
        # stdin whether it may clear a pre-existing output directory, and Snakemake
        # always creates the parent directory of a rule's outputs, so that directory
        # ALWAYS exists by the time stripenn runs; with no terminal attached the
        # prompt reads EOF and click aborts. But feeding it through a PIPE trades one
        # bug for another: `yes` is killed by SIGPIPE the moment stripenn stops
        # reading, and under `set -o pipefail` that kills the rule -- after a run that
        # had actually succeeded, whose outputs Snakemake then deleted as incomplete.
        # A redirect is not a pipeline, so pipefail never sees the writer.
        #
        # --norm {params.norm}: RAW counts, deliberately. Feeding stripenn the
        # ICE-balanced matrix returns zero stripes genome-wide -- it is an image
        # method, cooler's balancing leaves NaN in every unmappable bin (24.2% of the
        # chr1 matrix here), and NaN survives the clip so Canny finds no edges.
        # Measured on Naive_CTCF_rep1/chr1 at identical thresholds: balanced 0
        # stripes, raw 38. See the note in config.yaml.
        #
        # No `|| true`. Swallowing the exit status is what made an abort and a
        # genuine empty result look identical for eleven libraries.
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

        # Only reachable when stripenn exited 0 (set -e above); a crash can no longer
        # arrive here and be recorded as "no stripes". stripenn does write both
        # tables unconditionally, header included, even when the frame is empty, so
        # this is a guard rather than the normal path -- but the header must match
        # what it actually emits. The columns below are the post-`drop` set from
        # stripenn.py: total/num/start/end/x/y/h/w/medpixel are dropped and
        # Stripiness is appended, so the previous header here named nine columns that
        # do not exist in the file and omitted the one the summary sorts on.
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
