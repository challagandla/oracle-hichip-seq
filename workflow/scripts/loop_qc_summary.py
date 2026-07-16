"""
Aggregate per-sample QC into a single JSON + Markdown summary with categorical
flags. Replicate QC uses PASS / FAIL / DISCORDANT / NOT_ASSESSED instead of
pretending a single-replicate sample passed HiCRep or selecting its best match.
"""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from balance_utils import load_balance_report  # noqa: E402
from utils import load_loops_bedpe, setup_logging, write_json  # noqa: E402


DEFAULT_THRESHOLDS = {
    "valid_pair_yield_pct": 25.0,
    "duplicate_pct_max": 50.0,
    "cis_fraction_min": 0.70,
    # hicrep_scc_min is supplied from config at runtime (one source of truth).
}

HICREP_STATUSES = frozenset({"PASS", "FAIL", "DISCORDANT", "NOT_ASSESSED"})


def parse_pairtools_stats(path: str | Path) -> dict[str, float]:
    """Parse the key/value style file emitted by `pairtools stats`."""
    out: dict[str, float] = {}
    for line in Path(path).read_text().splitlines():
        m = re.match(r"^(\S+)\s+(\S+)", line)
        if not m:
            continue
        k, v = m.group(1), m.group(2)
        try:
            out[k] = float(v)
        except ValueError:
            pass
    return out


def _status(flag: bool) -> str:
    return "PASS" if flag else "FAIL"


def hicrep_status_decision(value: object) -> tuple[str, bool | None]:
    """Return the declared HiCRep category and its hard-gate interpretation.

    Numeric SCC summaries are deliberately absent from this function. They are
    descriptive values, while the upstream pairwise classifier is the sole
    source of the replicate-concordance decision.
    """
    status = str(value if value is not None else "NOT_ASSESSED")
    if status not in HICREP_STATUSES:
        raise ValueError(
            f"Unknown HiCRep status {status!r}; expected one of "
            f"{sorted(HICREP_STATUSES)}"
        )
    if status == "PASS":
        return status, True
    if status == "FAIL":
        return status, False
    return status, None


def parse_one_row_tsv(path: str | Path) -> dict[str, float | str]:
    lines = [line.split("\t") for line in Path(path).read_text().splitlines() if line]
    if len(lines) < 2:
        return {}
    row: dict[str, float | str] = {}
    for key, value in zip(lines[0], lines[1]):
        try:
            row[key] = float(value)
        except ValueError:
            row[key] = value
    return row


def fastp_pair_populations(path: str | Path) -> tuple[int, int]:
    """Return raw and retained read-pair counts from a paired-end fastp report."""
    report = json.loads(Path(path).read_text())
    values = []
    for stage in ("before_filtering", "after_filtering"):
        try:
            reads = int(report["summary"][stage]["total_reads"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                f"fastp JSON lacks integer summary.{stage}.total_reads: {path}"
            ) from exc
        if reads < 0 or reads % 2:
            raise ValueError(
                f"paired-end fastp total_reads must be a non-negative even number; "
                f"got {reads} at {stage}"
            )
        values.append(reads // 2)
    if values[1] > values[0]:
        raise ValueError("fastp retained more read pairs than it received")
    return values[0], values[1]


def contact_qc_metrics(
    pair_stats: dict[str, float],
    dedup_stats: dict[str, float],
    raw_input_pairs: int,
    post_trim_pairs: int,
) -> dict[str, float]:
    """Return metrics with explicit, non-overlapping denominators.

    ``pair_stats`` is calculated after selecting deduplicated, valid-ligation UU
    pairs, whereas ``dedup_stats`` begins after fastp and alignment. Fastp is the
    authoritative source for all sequenced read pairs and post-trim read pairs.
    Pairtools' own duplicate fraction is duplicates / mapped pairs; raw valid-pair
    yield is final high-confidence UU contacts / all sequenced read pairs.
    """
    pairtools_input_pairs = float(dedup_stats.get("total", 0.0))
    mapped_pairs = float(dedup_stats.get("total_mapped", 0.0))
    selected_uu = float(pair_stats.get("total", pair_stats.get("total_nodups", 0.0)))
    duplicate_pairs = float(dedup_stats.get("total_dups", 0.0))
    cis_pairs = float(pair_stats.get("cis", 0.0))

    duplicate_fraction = dedup_stats.get("summary/frac_dups")
    if duplicate_fraction is None:
        duplicate_fraction = (
            duplicate_pairs / mapped_pairs
            if mapped_pairs else 0.0
        )

    cis_fraction = pair_stats.get("summary/frac_cis")
    if cis_fraction is None:
        cis_fraction = (
            cis_pairs / selected_uu if selected_uu else 0.0
        )

    return {
        "raw_input_pairs": float(raw_input_pairs),
        "post_trim_pairs": float(post_trim_pairs),
        "pairtools_input_pairs": pairtools_input_pairs,
        "mapped_pairs": mapped_pairs,
        "duplicate_pairs": duplicate_pairs,
        "high_confidence_pairs": selected_uu,
        "cis_pairs": cis_pairs,
        "valid_pair_yield_pct": (
            100.0 * selected_uu / raw_input_pairs if raw_input_pairs else 0.0
        ),
        "post_trim_valid_pair_yield_pct": (
            100.0 * selected_uu / post_trim_pairs if post_trim_pairs else 0.0
        ),
        "trim_retention_pct": (
            100.0 * post_trim_pairs / raw_input_pairs if raw_input_pairs else 0.0
        ),
        "duplicate_pct": 100.0 * float(duplicate_fraction),
        "cis_fraction": float(cis_fraction),
    }


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])

    thresholds = dict(DEFAULT_THRESHOLDS)
    thresholds.update(dict(snakemake.params.thresholds))

    pair_stats = parse_pairtools_stats(snakemake.input.pair_stats)
    dedup_stats = parse_pairtools_stats(snakemake.input.dedup_stats)
    raw_input_pairs, post_trim_pairs = fastp_pair_populations(snakemake.input.fastp)
    apa = json.loads(Path(snakemake.input.apa_json).read_text())
    hicrep = json.loads(Path(snakemake.input.hicrep).read_text())
    balance = load_balance_report(snakemake.input.balance)
    contact_depth = json.loads(Path(snakemake.input.contact_depth).read_text())
    mustache_inputs = list(getattr(snakemake.input, "mustache", []))
    mustache = (
        json.loads(Path(mustache_inputs[0]).read_text())
        if mustache_inputs else None
    )
    restriction = json.loads(Path(snakemake.input.restriction).read_text())
    anchor_qc = parse_one_row_tsv(snakemake.input.anchor_qc)
    loops = load_loops_bedpe(snakemake.input.loops)

    metrics = contact_qc_metrics(
        pair_stats, dedup_stats, raw_input_pairs, post_trim_pairs
    )
    n_loops = len(loops)
    valid_yield = metrics["valid_pair_yield_pct"]
    dup_pct = metrics["duplicate_pct"]
    cis_frac = metrics["cis_fraction"]
    # APA is a descriptive held-out effect size. Even with anchor/visibility-
    # matched controls, its bootstrap interval is not a universally calibrated
    # assay-quality threshold and therefore never gates the sample.
    apa_matched_effect = apa.get("apa_vs_random_shift")
    apa_corner_ratio = apa.get("apa_score")
    apa_normalization = apa.get("normalization")
    # Numeric SCC summaries remain visible, but the categorical upstream decision
    # is authoritative. Re-thresholding the best match here would allow one strong
    # sibling pair to hide another depth-qualified discordant pair.
    hicrep_status, hicrep_pass = hicrep_status_decision(hicrep.get("status"))
    hicrep_group_status, hicrep_group_pass = hicrep_status_decision(
        hicrep.get("group_status")
    )
    hicrep_min = hicrep.get("min_scc")
    hicrep_mean = hicrep.get("mean_scc")
    hicrep_best = hicrep.get("best_scc")

    pass_flags = {
        "valid_pair_yield": valid_yield >= thresholds["valid_pair_yield_pct"],
        "duplicate_pct": dup_pct <= thresholds["duplicate_pct_max"],
        "cis_fraction": cis_frac >= thresholds["cis_fraction_min"],
    }
    status_flags = {k: _status(v) for k, v in pass_flags.items()}
    status_flags["n_loops"] = "DESCRIPTIVE"

    # Too few matched loops is NOT_ASSESSED rather than a false zero. Assessed APA
    # is DESCRIPTIVE and has no boolean pass/fail interpretation.
    apa_status = str(apa.get("status", "NOT_ASSESSED"))
    if apa_status not in {"DESCRIPTIVE", "NOT_ASSESSED"}:
        raise ValueError(f"Unknown APA status {apa_status!r}")
    status_flags["apa_matched_effect"] = apa_status
    # Kept as an alias for readers of the v1 report contract.
    status_flags["apa_score"] = apa_status
    status_flags["hicrep_scc"] = hicrep_status
    status_flags["hicrep_group"] = hicrep_group_status
    balance_status = str(balance["status"])
    status_flags["cooler_balance"] = balance_status
    mustache_status = "DISABLED"
    if mustache is not None:
        mustache_status = str(mustache.get("status", "NOT_ASSESSED"))
        if mustache_status not in {"PASS", "NOT_ASSESSED"}:
            raise ValueError(f"Unknown Mustache status {mustache_status!r}")
        status_flags["mustache"] = mustache_status

    legacy_frip = anchor_qc.get("frip")
    sample_peak_frip = anchor_qc.get("sample_peak_frip", legacy_frip)
    consensus_peak_frip = anchor_qc.get("consensus_peak_frip", legacy_frip)
    n_sample_peaks = int(
        anchor_qc.get("n_sample_peaks", anchor_qc.get("n_consensus_peaks", 0))
    )
    n_consensus_peaks = int(anchor_qc.get("n_consensus_peaks", 0))
    frip_population = str(
        anchor_qc.get("analysis_population", "primary_autosomes_chrX")
    )

    hard_pass = (
        all(pass_flags.values())
        and (hicrep_pass is not False)
        and (hicrep_group_pass is not False)
    )
    if not hard_pass:
        overall_status = "FAIL"
    elif {"DISCORDANT", "WARN"}.intersection(status_flags.values()):
        overall_status = "PASS_WITH_UNCERTAINTY"
    elif "NOT_ASSESSED" in status_flags.values():
        overall_status = "PASS_WITH_NOT_ASSESSED"
    else:
        overall_status = "PASS"

    report = {
        "sample": snakemake.wildcards.sample,
        # input_pairs is retained as a backwards-compatible alias, but now has
        # the scientifically correct raw-read-pair population.
        "input_pairs": int(metrics["raw_input_pairs"]),
        "raw_input_pairs": int(metrics["raw_input_pairs"]),
        "post_trim_pairs": int(metrics["post_trim_pairs"]),
        "pairtools_input_pairs": int(metrics["pairtools_input_pairs"]),
        "mapped_pairs": int(metrics["mapped_pairs"]),
        "duplicate_pairs": int(metrics["duplicate_pairs"]),
        "high_confidence_pairs": int(metrics["high_confidence_pairs"]),
        "cis_pairs": int(metrics["cis_pairs"]),
        "primary_cis_offdiagonal_contacts": int(
            contact_depth["primary_cis_offdiagonal_contacts"]
        ),
        "fithichip_distance_range_contacts": int(
            contact_depth["fithichip_distance_range_contacts"]
        ),
        "fithichip_distance_range_bp": contact_depth["fithichip_distance_range_bp"],
        "valid_pair_yield_pct": float(valid_yield),
        "post_trim_valid_pair_yield_pct": float(
            metrics["post_trim_valid_pair_yield_pct"]
        ),
        "trim_retention_pct": float(metrics["trim_retention_pct"]),
        "duplicate_pct": float(dup_pct),
        "cis_fraction": float(cis_frac),
        "n_loops": int(n_loops),
        # APA corner ratio is only approximately distance matched: the two 3x3
        # corner neighbourhoods span offsets within +/-2 bins.  The held-out
        # random-shift comparison is the preferred descriptive effect.
        "apa_score": apa_corner_ratio,
        "apa_corner_ratio_near_distance": apa_corner_ratio,
        "apa_matched_effect": apa_matched_effect,
        "apa_vs_random_shift": apa_matched_effect,
        "apa_vs_random_shift_ci95": apa.get("apa_vs_random_shift_ci95"),
        "apa_status": apa_status,
        "apa_normalization": apa_normalization,
        "apa_balance_status": apa.get("balance_status"),
        "hicrep_min_scc": hicrep_min,
        "hicrep_mean_scc": hicrep_mean,
        "hicrep_best_scc": hicrep_best,
        "hicrep_status": hicrep_status,
        "hicrep_group_median_scc": hicrep.get("group_median_scc"),
        "hicrep_group_status": hicrep_group_status,
        "balance_status": balance_status,
        "balance_qc": {
            "status": balance_status,
            "n_configured": balance.get("n_configured"),
            "n_converged": balance.get("n_converged"),
            "n_nonconverged": balance.get("n_nonconverged"),
            "n_missing": balance.get("n_missing"),
            "converged_resolutions_bp": balance.get("converged_resolutions_bp", []),
            "nonconverged_resolutions_bp": balance.get("nonconverged_resolutions_bp", []),
            "missing_resolutions_bp": balance.get("missing_resolutions_bp", []),
        },
        "mustache_status": mustache_status,
        "mustache_balance_status": (
            mustache.get("balance_status") if mustache is not None else None
        ),
        "mustache_caller_comparison": (
            {
                key: mustache.get(key)
                for key in (
                    "comparison_resolution_bp",
                    "comparison_tolerance_bins",
                    "n_primary_loops_on_grid",
                    "n_mustache_loops_on_grid",
                    "n_overlapping_loops_on_grid",
                    "caller_jaccard",
                    "primary_supported_by_mustache_fraction",
                    "mustache_supported_by_primary_fraction",
                    "comparison_is_gate",
                )
            }
            if mustache is not None else None
        ),
        "restriction_qc": {
            "population": restriction.get(
                "population", "post_dedup_pre_contact_filter_UU_pairs"
            ),
            "denominator_description": restriction.get(
                "denominator_description",
                "Deduplicated high-confidence UU pairs before restriction-artifact contact filtering",
            ),
            "total_deduplicated_uu_pairs": restriction.get(
                "total_deduplicated_uu_pairs"
            ),
            "fractions": restriction.get("fractions", {}),
        },
        "restriction_artifact_fractions": restriction.get("fractions", {}),
        "sample_peak_frip": sample_peak_frip,
        "consensus_peak_frip": consensus_peak_frip,
        "anchor_frip_population": frip_population,
        # v1 compatibility alias: anchor_frip always means consensus FRiP.
        "anchor_frip": consensus_peak_frip,
        "n_sample_peaks": n_sample_peaks,
        "n_consensus_peaks": n_consensus_peaks,
        "thresholds": thresholds,
        "pass_flags": pass_flags,
        "status_flags": status_flags,
        "overall_status": overall_status,
        # Uncertainty and missing assessment are explicit non-passes for automated
        # ingestion even when every assessed hard gate passed.
        "overall_pass": overall_status == "PASS",
    }
    write_json(report, snakemake.output.json)

    if apa_matched_effect is not None:
        apa_line = (
            "- Contact-map-held-out APA matched effect: "
            f"**{apa_matched_effect:.2f}** ({apa_status}; {apa_normalization}; "
            "descriptive, not a pass/fail gate)"
        )
    else:
        apa_line = (
            "- Contact-map-held-out APA: **NOT_ASSESSED** "
            "(too few sibling-donor loops with matched controls)"
        )
    if hicrep_best is None:
        hicrep_line = f"- HiCRep replicate concordance: **{hicrep_status}** (no depth-qualified pair)"
    else:
        hicrep_line = (
            f"- HiCRep replicate concordance: **{hicrep_status}** "
            f"(min/mean/best SCC {hicrep_min:.3f}/{hicrep_mean:.3f}/{hicrep_best:.3f}; "
            f"pair threshold ≥ {thresholds['hicrep_scc_min']})"
        )
    hicrep_group_median = hicrep.get("group_median_scc")
    hicrep_group_line = (
        f"- HiCRep replicate-group status: **{hicrep_group_status}** "
        + (
            f"(descriptive median SCC {hicrep_group_median:.3f})"
            if hicrep_group_median is not None
            else "(no depth-qualified pair)"
        )
    )
    sample_frip_text = (
        f"{sample_peak_frip:.3f}" if sample_peak_frip is not None else "NOT_ASSESSED"
    )
    consensus_frip_text = (
        f"{consensus_peak_frip:.3f}"
        if consensus_peak_frip is not None else "NOT_ASSESSED"
    )
    apa_corner_line = (
        "- APA near-distance corner ratio: "
        + (
            f"**{apa_corner_ratio:.2f}** (two corner neighbourhoods, offsets "
            "within +/-2 bins; descriptive)"
            if apa_corner_ratio is not None
            else "**NOT_ASSESSED**"
        )
    )
    md_lines = [
        f"# QC report — {snakemake.wildcards.sample}",
        "",
        f"- High-confidence UU contacts: **{int(metrics['high_confidence_pairs']):,}** / {int(metrics['raw_input_pairs']):,} raw input pairs",
        f"- Raw valid-pair yield: **{valid_yield:.1f}%** (threshold ≥ {thresholds['valid_pair_yield_pct']}%)",
        f"- Post-trim valid-pair yield: **{metrics['post_trim_valid_pair_yield_pct']:.1f}%** (descriptive; denominator {int(metrics['post_trim_pairs']):,} retained pairs)",
        f"- Duplicate % of mapped pairs: **{dup_pct:.1f}%** (threshold ≤ {thresholds['duplicate_pct_max']}%)",
        f"- Cis fraction of high-confidence pairs: **{cis_frac:.2f}** (threshold ≥ {thresholds['cis_fraction_min']:.2f})",
        f"- Significant FitHiChIP loops: **{n_loops}** (descriptive; not a universal QC gate)",
        f"- Read-end FRiP in this sample's MACS3 peaks: **{sample_frip_text}** ({n_sample_peaks:,} peaks; denominator {frip_population})",
        f"- Read-end FRiP in assay-stratum consensus peaks: **{consensus_frip_text}** ({n_consensus_peaks:,} peaks; denominator {frip_population})",
        "- Restriction-orientation denominator: **post-dedup high-confidence UU pairs before contact filtering**",
        apa_line,
        apa_corner_line,
        hicrep_line,
        hicrep_group_line,
        (
            f"- Cooler balancing: **{balance_status}** "
            f"({balance.get('n_converged', 0)}/{balance.get('n_configured', 0)} "
            "configured resolutions converged)"
        ),
        (
            f"- Mustache secondary cross-check: **{mustache_status}**"
            + (
                f" (balance {mustache.get('balance_status', 'NOT_ASSESSED')})"
                if mustache is not None else " (disabled by configuration)"
            )
        ),
        "",
        f"**Overall: {overall_status}**",
    ]
    if mustache is not None and mustache.get("n_mustache_loops_on_grid") is not None:
        md_lines.insert(
            -2,
            "- Caller overlap on the common "
            f"{mustache.get('comparison_resolution_bp', 0):,}-bp grid: "
            f"**{mustache.get('n_overlapping_loops_on_grid', 0):,}** shared / "
            f"{mustache.get('n_primary_loops_on_grid', 0):,} FitHiChIP / "
            f"{mustache.get('n_mustache_loops_on_grid', 0):,} Mustache "
            "(descriptive; not a gate)",
        )
    Path(snakemake.output.md).write_text("\n".join(md_lines))


# Guarded so the module can be imported by the tests. Snakemake injects
# `snakemake` into the script's globals before executing it.
if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
