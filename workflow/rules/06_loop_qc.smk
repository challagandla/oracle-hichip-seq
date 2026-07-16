# Stage 06 — Matrix and loop QC
# Cis/trans ratio, P(s) distance decay, exploratory local insulation and E1,
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

    chrY is excluded from the shared structural view to avoid ploidy- and
    mappability-dependent behavior and to keep the rule valid for replacement
    cohorts. chrX remains available for locus tracks and exploratory structural
    context; HiCRep separately uses an autosomal-only contract.
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
    P(s) distance-decay curve for checking a broadly decreasing contact
    frequency with genomic separation. HiChIP enrichment can alter its shape,
    so no universal slope is enforced.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        balance = RESULTS / "qc/balance/{sample}.balance.json",
        view = RESULTS / "qc/view_main_chroms.bed",
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        tsv = RESULTS / "qc/expected/{sample}.expected.cis.tsv"
    params:
        kind = "expected_cis",
        res = 25000,
        ignore_diags = config["cooler"]["balance"]["ignore_diags"],
    threads: 4
    conda: "../envs/cooltools.yaml"
    log:
        RESULTS / "logs/cooltools_expected/{sample}.log"
    script:
        "../scripts/cooltools_matrix_qc.py"


rule contact_depth_qc:
    """Depth denominators matched to loop and stripe caller search spaces."""
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        view = RESULTS / "qc/view_main_chroms.bed",
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        json = RESULTS / "qc/contact_depth/{sample}.json",
    params:
        resolution = config["fithichip"]["bin_size"],
        lower_distance = config["fithichip"]["lower_distance"],
        upper_distance = config["fithichip"]["upper_distance"],
    conda: "../envs/coolerpy.yaml"
    log:
        RESULTS / "logs/contact_depth/{sample}.log",
    script:
        "../scripts/contact_depth_qc.py"

rule cooltools_insulation:
    """Exploratory local-insulation score at 25 kb.

    HiChIP contact maps are enriched around the immunoprecipitated mark, so this
    is useful for within-mark locus context but is not an unbiased genome-wide
    TAD-boundary call. Canonical boundary claims require orthogonal Hi-C/Micro-C
    evidence or a design that explicitly validates the HiChIP-derived signal.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        balance = RESULTS / "qc/balance/{sample}.balance.json",
        view = RESULTS / "qc/view_main_chroms.bed",
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        tsv = RESULTS / "qc/insulation/{sample}.insulation.tsv"
    params:
        kind = "insulation",
        res = 25000,
        window = 250000,
        ignore_diags = config["cooler"]["balance"]["ignore_diags"],
    threads: 4
    conda: "../envs/cooltools.yaml"
    log:
        RESULTS / "logs/cooltools_insulation/{sample}.log"
    script:
        "../scripts/cooltools_matrix_qc.py"

rule gc_phasing_track:
    """GC content per 100 kb bin, used to orient the exploratory E1 vector.

    The sign of an eigenvector is arbitrary. cooltools solves E1 per chromosome, so
    without a phasing track "A" is positive on some chromosomes and negative on
    others, independently in every sample. Correlating raw E1 between two libraries
    can therefore average arbitrary sign flips towards zero. This is an orientation
    problem; fixing it does not by itself establish biological compartment validity.

    GC content is a standard phasing track -- the A compartment is generally the
    gene-rich, GC-rich one -- so it fixes the sign consistently. Phasing resolves
    the arbitrary sign only; it does not make a mark-enriched HiChIP matrix an
    unbiased compartment assay or remove depth and antibody-enrichment effects.
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
    """GC-phased E1-like signal at 100 kb, normalised to a stable TSV schema.

    Interpret within the same mark and adequate-depth libraries as exploratory
    domain-scale context. Do not call canonical A/B compartment differences from
    HiChIP alone; confirm them with Hi-C/Micro-C or another unbiased contact assay.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        balance = RESULTS / "qc/balance/{sample}.balance.json",
        view = RESULTS / "qc/view_main_chroms.bed",
        gc = RESULTS / "qc/compartments/gc_100kb.tsv",
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        cis = RESULTS / "qc/compartments/{sample}.cis.eigs.tsv",
        status = RESULTS / "qc/compartments/{sample}.cis.eigs.status.json",
    params:
        res = 100000
    threads: 4
    conda: "../envs/cooltools.yaml"
    log:
        RESULTS / "logs/cooltools_eigs/{sample}.log"
    script:
        "../scripts/cooltools_eigs_cis.py"

rule compartments_to_bigwig:
    """Export exploratory GC-phased E1 to a browser/pyGenomeTracks bigWig."""
    input:
        eigs = RESULTS / "qc/compartments/{sample}.cis.eigs.tsv",
        chromsizes = GENOME["chromsizes"],
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        bw = RESULTS / "qc/compartments/{sample}.E1.bw"
    threads: 1
    conda: "../envs/coolerpy.yaml"
    log:
        RESULTS / "logs/compartments_to_bigwig/{sample}.log"
    script:
        "../scripts/compartments_to_bigwig.py"

HICREP_GROUP_BY = list(config["hicrep"]["group_by"])


def _hicrep_group_samples(sample):
    row = SAMPLES.loc[sample]
    selected = SAMPLES.copy()
    for column in HICREP_GROUP_BY:
        selected = selected[
            selected[column].fillna("").astype(str)
            == str(row[column] if pd.notna(row[column]) else "")
        ]
    return selected["sample_id"].tolist()


rule hicrep_replicate_qc:
    """
    Stratum-adjusted correlation between biological replicates in the same
    configured condition/mark/protocol group and different donors. Depth-qualified
    pairwise results are classified as all-pass, all-fail, or discordant; no
    best-pair shortcut is used.

    Replicates are grouped by the configured biological condition, mark, tissue,
    and protocol columns, never by subject_id. Grouping on donor would put distinct
    conditions from the same person into one "replicate" group and score biological
    differences against a threshold intended for technical concordance.
    """
    input:
        mcools = lambda wc: expand(
            RESULTS / "cool/{sample}.mcool",
            sample=_hicrep_group_samples(wc.sample),
        ),
        view = RESULTS / "qc/view_main_chroms.bed",
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        json = RESULTS / "qc/hicrep/{sample}.hicrep.json"
    params:
        bin = config["hicrep"]["bin_size"],
        maxd = config["hicrep"]["max_dist"],
        h = config["hicrep"]["h_smooth"],
        threshold = config["hicrep"]["threshold_pass"],
        min_contacts = config["hicrep"]["min_contacts_for_scc"],
    threads: 4
    conda: "../envs/hicrep.yaml"
    log:
        RESULTS / "logs/hicrep/{sample}.log"
    script:
        "../scripts/hicrep_replicate_qc.py"

def _sibling_sample_ids(wc):
    row = SAMPLES.loc[wc.sample]
    replicate_group = SAMPLES.loc[_hicrep_group_samples(wc.sample)]
    return replicate_group[
        (replicate_group["subject_id"] != row["subject_id"]) &
        (replicate_group["sample_id"] != wc.sample)
    ]["sample_id"].tolist()


def _heldout_loop_inputs(wc):
    return expand(
        RESULTS / f"loops/{{sample}}/{{sample}}.interactions_FitHiC_{FITHICHIP_Q_LABEL}.bed",
        sample=_sibling_sample_ids(wc),
    )


def _heldout_peak_inputs(wc):
    return expand(
        RESULTS / "peaks/{sample}_peaks.bed",
        sample=_sibling_sample_ids(wc),
    )


def _heldout_support(wc):
    # A three-replicate group has two held-out donors and requires agreement between
    # both; a two-replicate group has one held-out call set and can require only one.
    return min(2, max(1, len(_heldout_loop_inputs(wc))))


rule heldout_apa_loops:
    """Build sibling-donor candidates without using the map or peaks being scored.

    The sibling FitHiChIP calls were generated in the shared configured assay-stratum anchor
    universe. Requiring overlap with a sibling-only consensus anchor removes
    target-only anchor leakage. The residual shared search-space dependence is
    retained explicitly in the APA metadata rather than described as strict
    statistical independence.
    """
    input:
        loops = _heldout_loop_inputs,
        blacklist = GENOME["blacklist"],
        chromsizes = GENOME["chromsizes"],
        anchors = RESULTS / "qc/apa_candidates/{sample}.heldout_anchors.bed",
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        bedpe = RESULTS / "qc/apa_candidates/{sample}.heldout.bedpe",
        audit = RESULTS / "qc/apa_candidates/{sample}.heldout_support.tsv",
    params:
        bin_size = config["fithichip"]["bin_size"],
        min_sample_support = _heldout_support,
        tolerance_bins = config["apa"]["candidate_tolerance_bins"],
        allow_empty = True,
    conda: "../envs/pandas.yaml"
    log:
        RESULTS / "logs/apa_candidates/{sample}.log",
    script:
        "../scripts/build_consensus_loops.py"


rule heldout_apa_peaks:
    """Build a same-group anchor set from sibling donors only."""
    input:
        peaks = _heldout_peak_inputs,
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        bed = RESULTS / "qc/apa_candidates/{sample}.heldout_anchors.bed",
        audit = RESULTS / "qc/apa_candidates/{sample}.heldout_anchor_support.tsv",
    params:
        min_support = _heldout_support,
        allow_empty = True,
    conda: "../envs/pandas.yaml"
    log:
        RESULTS / "logs/apa_candidates/{sample}.anchors.log",
    script:
        "../scripts/consensus_peaks.py"


rule apa_plot:
    """
    Contact-map-held-out Aggregate Peak Analysis: score this map on loops reproduced
    by sibling donors and supported by sibling-only anchors, with distance-matched
    random-shift controls and uncertainty.
    """
    input:
        mcool = RESULTS / "cool/{sample}.mcool",
        balance = RESULTS / "qc/balance/{sample}.balance.json",
        loops = RESULTS / "qc/apa_candidates/{sample}.heldout.bedpe",
        candidate_anchors = RESULTS / "qc/apa_candidates/{sample}.heldout_anchors.bed",
        candidate_audit = RESULTS / "qc/apa_candidates/{sample}.heldout_support.tsv",
        candidate_anchor_audit = RESULTS / "qc/apa_candidates/{sample}.heldout_anchor_support.tsv",
        blacklist = GENOME["blacklist"],
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        png = RESULTS / "qc/apa/{sample}.apa.png",
        json = RESULTS / "qc/apa/{sample}.apa.json",
        npy = RESULTS / "qc/apa/{sample}.apa.npy"
    params:
        window = config["apa"]["window_size"],
        bin_size = config["apa"]["bin_size"],
        min_dist = config["apa"]["min_loop_dist"],
        n_ctrl = config["apa"]["n_random_controls"],
        min_loops = config["apa"]["min_loops_for_apa"],
        control_marginal_log2_tolerance = config["apa"]["control_marginal_log2_tolerance"],
        max_control_attempts_per_draw = config["apa"]["max_control_attempts_per_draw"],
        visibility_min_dist = config["fithichip"]["lower_distance"],
        visibility_max_dist = config["fithichip"]["upper_distance"],
        n_sibling_donor_callsets = lambda wc: len(_heldout_loop_inputs(wc)),
        candidate_min_sample_support = _heldout_support,
        candidate_tolerance_bins = config["apa"]["candidate_tolerance_bins"],
        candidate_grid_bin_size_bp = config["fithichip"]["bin_size"],
        contact_map_held_out = True,
        anchors_exclude_scored_sample = True,
        primary_call_search_space = (
            "shared configured assay-stratum consensus anchors include the scored sample; "
            "q-filtered sibling calls are reconciled within the configured APA bin "
            "tolerance and subsequently restricted to sibling-only anchors; this is "
            "separate from the exact differential hypothesis grid"
        ),
    threads: 4
    conda: "../envs/coolerpy.yaml"
    log:
        RESULTS / "logs/apa/{sample}.log"
    script:
        "../scripts/apa_plot.py"

rule loop_qc_summary:
    """
    Aggregate every QC metric into one JSON per sample with explicit pass, fail,
    discordant, descriptive, and not-assessed states.
    Consumed by MultiQC custom content.
    """
    input:
        pair_stats = RESULTS / "qc/pairtools/{sample}.pairs.stats.txt",
        dedup_stats = RESULTS / "qc/pairtools/{sample}.dedup.stats.txt",
        fastp = RESULTS / "qc/fastp/{sample}.fastp.json",
        restriction = RESULTS / "qc/restriction/{sample}.restriction.json",
        anchor_qc = RESULTS / "qc/anchors/{sample}.anchor_qc.tsv",
        expected = RESULTS / "qc/expected/{sample}.expected.cis.tsv",
        apa_json = RESULTS / "qc/apa/{sample}.apa.json",
        hicrep = RESULTS / "qc/hicrep/{sample}.hicrep.json",
        balance = RESULTS / "qc/balance/{sample}.balance.json",
        contact_depth = RESULTS / "qc/contact_depth/{sample}.json",
        mustache = lambda wc: (
            [RESULTS / f"loops/{wc.sample}/{wc.sample}.mustache.status.json"]
            if config.get("mustache", {}).get("enabled") else []
        ),
        loops = RESULTS / f"loops/{{sample}}/{{sample}}.interactions_FitHiC_{FITHICHIP_Q_LABEL}.bed",
        shared_code = SHARED_SCRIPT_DEPS,
    output:
        json = RESULTS / "qc/loop_qc/{sample}.json",
        md   = RESULTS / "qc/loop_qc/{sample}.md"
    threads: 1
    params:
        thresholds = {
            **config.get("qc_thresholds", {}),
            "hicrep_scc_min": config["hicrep"]["threshold_pass"],
        },
    conda: "../envs/pandas.yaml"
    log:
        RESULTS / "logs/loop_qc_summary/{sample}.log"
    script:
        "../scripts/loop_qc_summary.py"
