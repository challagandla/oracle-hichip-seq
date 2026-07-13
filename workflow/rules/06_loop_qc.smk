# Stage 06 — Matrix and loop QC
# Cis/trans ratio, P(s) distance decay, insulation scores, A/B compartments,
# APA aggregate analysis on the called loop set, HiCRep replicate concordance.


rule main_chrom_view:
    """A cooltools 'view' restricted to the assembled chromosomes.

    hg38 carries ~160 unplaced scaffolds and alt contigs. cooltools and HiCRep
    otherwise iterate every region in the cooler, and on a scaffold with no valid
    bins after balancing they do not skip it -- they die:

        cooltools insulation  IndexError: index 0 is out of bounds for axis 0 with size 0
        hicrepSCC             AssertionError: Contact matrix 1 of chromosome GL000208.1 is empty

    Restricting them is not just a workaround. Insulation, compartments and
    stratum-adjusted correlation are all defined on a chromosome with a real
    distance-decay profile; on a 60 kb unplaced contig they are meaningless even
    when they happen to compute.

    chrY is dropped as well: donor sex is not recorded in this cohort, so it is
    present in some libraries and absent in others, and a region that exists for
    only some samples cannot be compared across them.
    """
    input:
        chromsizes = GENOME["chromsizes"],
    output:
        view = RESULTS / "qc/view_main_chroms.bed",
    conda: "../envs/coreutils.yaml"
    log:
        RESULTS / "logs/main_chrom_view.log",
    shell:
        r"""
        mkdir -p $(dirname {output.view}) $(dirname {log})
        awk 'BEGIN{{OFS="\t"}} $1 ~ /^chr([0-9]+|X)$/ {{print $1, 0, $2, $1}}' \
            {input.chromsizes} | sort -k1,1V > {output.view} 2> {log}
        test -s {output.view}
        echo "view regions: $(wc -l < {output.view})" >> {log}
        """


rule cooltools_expected_cis:
    """
    P(s) distance-decay curve. The expected −1 slope on log–log is the
    canonical sanity check for any HiC-style assay.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        view = RESULTS / "qc/view_main_chroms.bed",
    output:
        expected = RESULTS / "qc/expected/{sample}.expected.cis.tsv"
    params:
        res = 25000
    threads: 4
    conda: "../envs/cooltools.yaml"
    log:
        RESULTS / "logs/cooltools_expected/{sample}.log"
    shell:
        r"""
        cooltools expected-cis -p {threads} --view {input.view} \
            {input.mcool}::resolutions/{params.res} \
            > {output.expected} 2> {log}
        """

rule cooltools_insulation:
    """Insulation score (TAD boundary detection)."""
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        view = RESULTS / "qc/view_main_chroms.bed",
    output:
        tsv = RESULTS / "qc/insulation/{sample}.insulation.tsv"
    params:
        res = 25000,
        window = 250000
    threads: 4
    conda: "../envs/cooltools.yaml"
    log:
        RESULTS / "logs/cooltools_insulation/{sample}.log"
    shell:
        r"""
        # --view: without it cooltools walks every unplaced scaffold in the cooler
        # and dies on the first one with no valid bins after balancing
        # (IndexError: index 0 is out of bounds for axis 0 with size 0).
        cooltools insulation \
            -p {threads} \
            --view {input.view} \
            {input.mcool}::resolutions/{params.res} \
            {params.window} \
            > {output.tsv} 2> {log}
        """

rule gc_phasing_track:
    """GC content per 100 kb bin, used to orient the compartment eigenvector.

    The sign of an eigenvector is arbitrary. cooltools solves E1 per chromosome, so
    without a phasing track "A" is positive on some chromosomes and negative on
    others, independently in every sample. Correlating raw E1 between two libraries
    then averages random sign flips towards zero: two Naive replicates that HiCRep
    scores at 0.867 came out at 0.20 on compartment agreement, and some genuine
    replicate pairs were NEGATIVE. The compartments were right; only their
    orientation was unfixed.

    GC content is the standard phasing track -- the A compartment is the
    gene-rich, GC-rich one -- so it fixes the sign the same way in every sample and
    makes E1 comparable across libraries.
    """
    input:
        chromsizes = GENOME["chromsizes"],
        fasta = GENOME["fasta"],
        view = RESULTS / "qc/view_main_chroms.bed",
    output:
        gc = RESULTS / "qc/compartments/gc_100kb.tsv",
    params:
        res = 100000,
    threads: 1
    conda: "../envs/cooltools.yaml"
    log:
        RESULTS / "logs/gc_phasing_track.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.gc}) $(dirname {log})
        tmp=$(dirname {output.gc})/.bins_{params.res}.bed
        # Binnify the assembled chromosomes only, matching the view every other
        # cooltools rule uses; a GC value for an unplaced scaffold has nothing to
        # phase against.
        cooltools genome binnify --all-names {input.chromsizes} {params.res} \
            | awk 'NR==1 || $1 ~ /^chr([0-9]+|X)$/' > "$tmp" 2> {log}
        cooltools genome gc "$tmp" {input.fasta} > {output.gc} 2>> {log}
        rm -f "$tmp"
        test -s {output.gc}
        """


rule cooltools_eigs_cis:
    """A/B compartment eigenvectors at 100 kb, normalised to a stable TSV schema."""
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        view = RESULTS / "qc/view_main_chroms.bed",
        gc = RESULTS / "qc/compartments/gc_100kb.tsv",
    output:
        cis = RESULTS / "qc/compartments/{sample}.cis.eigs.tsv"
    params:
        res = 100000
    threads: 4
    conda: "../envs/cooltools.yaml"
    log:
        RESULTS / "logs/cooltools_eigs/{sample}.log"
    script:
        "../scripts/cooltools_eigs_cis.py"

rule compartments_to_bigwig:
    """Export E1 (A/B compartments) to a bigWig for browser / pyGenomeTracks."""
    input:
        eigs = RESULTS / "qc/compartments/{sample}.cis.eigs.tsv"
    output:
        bw = RESULTS / "qc/compartments/{sample}.E1.bw"
    params:
        chromsizes = GENOME["chromsizes"]
    threads: 1
    conda: "../envs/coolerpy.yaml"
    log:
        RESULTS / "logs/compartments_to_bigwig/{sample}.log"
    script:
        "../scripts/compartments_to_bigwig.py"

rule hicrep_replicate_qc:
    """
    Stratum-adjusted correlation between biological replicates: same cell type,
    same mark, different donors. The only concordance metric robust to
    distance-decay differences between libraries.

    Replicates are grouped on cell_type + mark, NOT subject_id + mark. Grouping on
    the donor puts Naive, Th17 and Treg from the same person into one "replicate"
    group and reports the correlation BETWEEN CELL TYPES as replicate concordance
    -- which is a measurement of the biology the differential test is trying to
    find, scored against a threshold that assumes it is measuring noise.
    """
    input:
        mcools = lambda wc: expand(
            RESULTS / "cool/{sample}.mcool",
            sample=SAMPLES[
                (SAMPLES["cell_type"] == SAMPLES.loc[wc.sample, "cell_type"]) &
                (SAMPLES["mark"] == SAMPLES.loc[wc.sample, "mark"])
            ]["sample_id"].tolist()
        ),
        view = RESULTS / "qc/view_main_chroms.bed"
    output:
        json = RESULTS / "qc/hicrep/{sample}.hicrep.json"
    params:
        bin = config["hicrep"]["bin_size"],
        maxd = config["hicrep"]["max_dist"],
        h = config["hicrep"]["h_smooth"]
    threads: 4
    conda: "../envs/hicrep.yaml"
    log:
        RESULTS / "logs/hicrep/{sample}.log"
    script:
        "../scripts/hicrep_replicate_qc.py"

rule apa_plot:
    """
    Aggregate Peak Analysis on the called loop set with random-shift controls.
    A score ≥ 1.5 vs controls is the standard quality cutoff.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        loops = RESULTS / f"loops/{{sample}}/{{sample}}.interactions_FitHiC_{FITHICHIP_Q_LABEL}.bed"
    output:
        png = RESULTS / "qc/apa/{sample}.apa.png",
        json = RESULTS / "qc/apa/{sample}.apa.json",
        npy = RESULTS / "qc/apa/{sample}.apa.npy"
    params:
        window = config["apa"]["window_size"],
        bin_size = config["apa"]["bin_size"],
        min_dist = config["apa"]["min_loop_dist"],
        n_ctrl = config["apa"]["n_random_controls"]
    threads: 4
    conda: "../envs/coolerpy.yaml"
    log:
        RESULTS / "logs/apa/{sample}.log"
    script:
        "../scripts/apa_plot.py"

rule loop_qc_summary:
    """
    Aggregate every QC metric into a single JSON per sample with pass/fail/not-assessed flags.
    Consumed by MultiQC custom content.
    """
    input:
        pair_stats = RESULTS / "qc/pairtools/{sample}.pairs.stats.txt",
        dedup_stats = RESULTS / "qc/pairtools/{sample}.dedup.stats.txt",
        expected = RESULTS / "qc/expected/{sample}.expected.cis.tsv",
        apa_json = RESULTS / "qc/apa/{sample}.apa.json",
        hicrep = RESULTS / "qc/hicrep/{sample}.hicrep.json",
        loops = RESULTS / f"loops/{{sample}}/{{sample}}.interactions_FitHiC_{FITHICHIP_Q_LABEL}.bed"
    output:
        json = RESULTS / "qc/loop_qc/{sample}.json",
        md   = RESULTS / "qc/loop_qc/{sample}.md"
    threads: 1
    conda: "../envs/pandas.yaml"
    log:
        RESULTS / "logs/loop_qc_summary/{sample}.log"
    script:
        "../scripts/loop_qc_summary.py"
