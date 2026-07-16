"""Create self-contained MultiQC custom-content files from stable QC JSON.

Keeping this presentation step separate prevents a report-schema change from
re-running expensive APA calculations. The input JSON files remain the stable
machine-readable workflow contract.
"""
import json
from pathlib import Path


def apa_payload(report: dict) -> dict:
    """Build a MultiQC bar plot for one assessed held-out APA result."""
    score = report.get("apa_vs_random_shift")
    data = {}
    if score is not None and report.get("normalization") == "ICE-balanced":
        data[report["sample"]] = {"APA versus random shift": score}
    return {
        "id": "apa_scores",
        "section_name": "Held-out APA matched effects",
        "description": (
            "Descriptive contact-map-held-out Aggregate Peak Analysis effect at "
            "sibling-donor loop centres versus chromosome/distance, sibling-anchor-"
            "class, visibility, blacklist, and mappability-matched random shifts. "
            "Raw-count fallback values remain in APA JSON but are not mixed into "
            "this balanced-only plot."
        ),
        "plot_type": "bargraph",
        "pconfig": {
            "id": "apa_bar",
            "title": "Held-out APA centre / matched-control centre",
            "ylab": "Matched effect (fold)",
            "cpswitch": False,
        },
        # NOT_ASSESSED must not become a false numeric zero.
        "data": data,
    }


def balance_payload(report: dict) -> dict:
    """Build a self-contained table of Cooler convergence attributes."""
    fields = {
        "status": report.get("status", "NOT_ASSESSED"),
        "n_configured": report.get("n_configured", 0),
        "n_converged": report.get("n_converged", 0),
        "n_nonconverged": report.get("n_nonconverged", 0),
        "n_missing": report.get("n_missing", 0),
        "converged_resolutions": ", ".join(
            str(value) for value in report.get("converged_resolutions_bp", [])
        ) or "none",
        "nonconverged_resolutions": ", ".join(
            str(value) for value in report.get("nonconverged_resolutions_bp", [])
        ) or "none",
        "missing_resolutions": ", ".join(
            str(value) for value in report.get("missing_resolutions_bp", [])
        ) or "none",
    }
    return {
        "id": "balance_qc",
        "section_name": "Cooler Balance QC",
        "description": (
            "ICE convergence read directly from each bins/weight HDF5 dataset. "
            "Nonconverged weights are removed from the published mcool; supported "
            "downstream analyses use explicitly labelled raw-count fallbacks."
        ),
        "plot_type": "table",
        "pconfig": {"id": "balance_qc_table", "title": "Cooler balance convergence"},
        "headers": {
            "status": {"title": "Balance status"},
            "n_configured": {"title": "Configured", "format": "{:,.0f}"},
            "n_converged": {"title": "Converged", "format": "{:,.0f}"},
            "n_nonconverged": {"title": "Nonconverged", "format": "{:,.0f}"},
            "n_missing": {"title": "Missing", "format": "{:,.0f}"},
            "converged_resolutions": {"title": "Converged bp"},
            "nonconverged_resolutions": {"title": "Nonconverged bp"},
            "missing_resolutions": {"title": "Missing bp"},
        },
        "data": {report["sample"]: fields},
    }


def loop_qc_payload(report: dict) -> dict:
    """Build a MultiQC custom-content table for one loop-QC report."""
    thresholds = report["thresholds"]
    mustache_comparison = report.get("mustache_caller_comparison") or {}
    fields = {
        "valid_pair_yield_pct": report["valid_pair_yield_pct"],
        "post_trim_valid_pair_yield_pct": report.get(
            "post_trim_valid_pair_yield_pct"
        ),
        "duplicate_pct": report["duplicate_pct"],
        "cis_fraction": report["cis_fraction"],
        "sample_peak_frip": report.get(
            "sample_peak_frip", report.get("anchor_frip")
        ),
        "consensus_peak_frip": report.get(
            "consensus_peak_frip", report.get("anchor_frip")
        ),
        # Compatibility alias for existing MultiQC consumers.
        "anchor_frip": report.get(
            "consensus_peak_frip", report.get("anchor_frip")
        ),
        "n_sample_peaks": report.get(
            "n_sample_peaks", report.get("n_consensus_peaks")
        ),
        "n_consensus_peaks": report.get("n_consensus_peaks"),
        "frip_population": report.get(
            "anchor_frip_population", "primary_autosomes_chrX"
        ),
        "n_loops": report["n_loops"],
        "apa_matched_effect": report.get(
            "apa_matched_effect",
            report.get("apa_vs_random_shift", report.get("apa_score")),
        ),
        "apa_corner_ratio": report.get(
            "apa_corner_ratio_near_distance", report.get("apa_score")
        ),
        "apa_status": report.get("apa_status", "NOT_ASSESSED"),
        "apa_normalization": report.get("apa_normalization"),
        "hicrep_status": report.get("hicrep_status", "NOT_ASSESSED"),
        "hicrep_min_scc": report.get("hicrep_min_scc"),
        "hicrep_mean_scc": report.get("hicrep_mean_scc"),
        "hicrep_best_scc": report.get("hicrep_best_scc"),
        "hicrep_group_status": report.get("hicrep_group_status", "NOT_ASSESSED"),
        "hicrep_group_median_scc": report.get("hicrep_group_median_scc"),
        "balance_status": report.get("balance_status", "NOT_ASSESSED"),
        "balance_converged": report.get("balance_qc", {}).get("n_converged"),
        "balance_configured": report.get("balance_qc", {}).get("n_configured"),
        "mustache_status": report.get("mustache_status", "DISABLED"),
        "mustache_balance_status": report.get("mustache_balance_status"),
        "mustache_n_loops": mustache_comparison.get("n_mustache_loops_on_grid"),
        "mustache_overlap": mustache_comparison.get("n_overlapping_loops_on_grid"),
        "mustache_jaccard": mustache_comparison.get("caller_jaccard"),
        "mustache_primary_supported": mustache_comparison.get(
            "primary_supported_by_mustache_fraction"
        ),
        "restriction_population": (report.get("restriction_qc") or {}).get(
            "population", "post_dedup_pre_contact_filter_UU_pairs"
        ),
        "overall_status": report["overall_status"],
    }
    return {
        "id": "loop_qc",
        "section_name": "Loop QC Summary",
        "description": (
            "Per-sample HiChIP QC covering contacts, assay-stratum anchors, "
            "loops, held-out APA, and replicate concordance."
        ),
        "plot_type": "table",
        "pconfig": {"id": "loop_qc_table", "title": "Loop QC"},
        "headers": {
            "valid_pair_yield_pct": {
                "title": "Raw valid-pair yield (%)",
                "description": "Final deduplicated valid-ligation UU pairs / raw sequenced read pairs x 100",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "post_trim_valid_pair_yield_pct": {
                "title": "Post-trim valid-pair yield (%)",
                "description": "Final deduplicated valid-ligation UU pairs / fastp-retained read pairs x 100; descriptive",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Blues",
            },
            "duplicate_pct": {
                "title": "Duplicate (%)",
                "description": "Duplicate pairs / mapped pairs x 100",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdYlGn-rev",
            },
            "cis_fraction": {
                "title": "Cis fraction",
                "description": "Cis high-confidence pairs / all high-confidence pairs",
                "min": 0,
                "max": 1,
                "scale": "RdYlGn",
            },
            "sample_peak_frip": {
                "title": "Sample-peak FRiP",
                "description": "Deduplicated UU read-end fraction in this library's own MACS3 peaks",
                "min": 0,
                "max": 1,
                "scale": "RdYlGn",
            },
            "consensus_peak_frip": {
                "title": "Consensus-peak FRiP",
                "description": "Deduplicated UU read-end fraction in assay-stratum consensus peaks",
                "min": 0,
                "max": 1,
                "scale": "RdYlGn",
            },
            "anchor_frip": {
                "title": "Consensus FRiP (v1 alias)",
                "description": "Compatibility alias of consensus-peak FRiP",
                "min": 0,
                "max": 1,
                "scale": "RdYlGn",
                "hidden": True,
            },
            "n_sample_peaks": {
                "title": "Sample peaks",
                "description": "Primary-contig, blacklist-filtered MACS3 peaks in this library",
                "format": "{:,.0f}",
                "scale": "Blues",
            },
            "n_consensus_peaks": {
                "title": "Consensus peaks",
                "description": "Cohort-supported anchor peaks",
                "format": "{:,.0f}",
                "scale": "Blues",
            },
            "frip_population": {
                "title": "FRiP population",
                "description": "Numerator and denominator both use assembled primary autosomes plus chrX",
            },
            "n_loops": {
                "title": "Significant loops",
                "description": "Filtered FitHiChIP loops at the configured FDR",
                "format": "{:,.0f}",
                "scale": "Blues",
            },
            "apa_matched_effect": {
                "title": "Held-out APA matched effect",
                "description": (
                    "Centre signal divided by matched random-shift centre signal; "
                    "descriptive and never a hard "
                    "sample pass/fail gate"
                ),
                "min": 0,
                "scale": "RdYlGn",
            },
            "apa_corner_ratio": {
                "title": "APA near-distance corner ratio",
                "description": (
                    "Centre divided by two 3x3 corner neighbourhoods whose "
                    "separation offsets are within +/-2 bins; descriptive"
                ),
                "min": 0,
                "scale": "Blues",
            },
            "apa_status": {
                "title": "APA status",
                "description": "DESCRIPTIVE when assessed; NOT_ASSESSED when matched evidence is insufficient",
            },
            "apa_normalization": {
                "title": "APA normalization",
                "description": "ICE-balanced or explicitly labelled raw-count fallback",
            },
            "hicrep_status": {
                "title": "HiCRep status",
                "description": (
                    "PASS when every depth-qualified pair involving this sample "
                    f"has SCC >= {thresholds['hicrep_scc_min']}; FAIL when every "
                    "pair is below; DISCORDANT when results are mixed"
                ),
            },
            "hicrep_min_scc": {
                "title": "HiCRep min SCC",
                "description": "Descriptive minimum across depth-qualified replicate pairs; not a gate",
                "min": -1,
                "max": 1,
                "scale": "RdYlGn",
            },
            "hicrep_mean_scc": {
                "title": "HiCRep mean SCC",
                "description": "Descriptive mean across depth-qualified replicate pairs; not a gate",
                "min": -1,
                "max": 1,
                "scale": "RdYlGn",
            },
            "hicrep_best_scc": {
                "title": "HiCRep best SCC",
                "description": "Descriptive maximum across depth-qualified replicate pairs; not a gate",
                "min": -1,
                "max": 1,
                "scale": "RdYlGn",
            },
            "hicrep_group_status": {
                "title": "HiCRep group",
                "description": (
                    "PASS when every depth-qualified pair in the replicate group "
                    "passes; FAIL when every pair fails; DISCORDANT when mixed"
                ),
            },
            "hicrep_group_median_scc": {
                "title": "HiCRep group median",
                "description": "Descriptive median across the replicate group; not a gate",
                "min": -1,
                "max": 1,
                "scale": "RdYlGn",
            },
            "balance_status": {
                "title": "Balance status",
                "description": "PASS when every configured resolution converged; WARN for any nonconverged weight",
            },
            "balance_converged": {
                "title": "Balanced resolutions",
                "description": "Configured resolutions whose bins/weight converged attribute is true",
                "format": "{:,.0f}",
            },
            "balance_configured": {
                "title": "Configured resolutions",
                "format": "{:,.0f}",
            },
            "mustache_status": {
                "title": "Mustache status",
                "description": "PASS when the balanced secondary cross-check ran; NOT_ASSESSED when its required balance was unavailable",
            },
            "mustache_balance_status": {
                "title": "Mustache balance",
                "description": "Convergence status at the configured Mustache resolution",
            },
            "mustache_n_loops": {
                "title": "Mustache loops",
                "description": "Secondary-caller loops projected to the common comparison grid",
                "format": "{:,.0f}",
            },
            "mustache_overlap": {
                "title": "Shared loops",
                "description": "One-to-one FitHiChIP/Mustache overlap within the configured reciprocal-anchor grid tolerance; descriptive, not a gate",
                "format": "{:,.0f}",
            },
            "mustache_jaccard": {
                "title": "Caller Jaccard",
                "description": "One-to-one tolerant overlap divided by the caller union; descriptive, not a gate",
                "min": 0,
                "max": 1,
                "scale": "Blues",
            },
            "mustache_primary_supported": {
                "title": "FitHiChIP supported",
                "description": "Fraction of FitHiChIP grid loops also called by Mustache; descriptive, not a gate",
                "min": 0,
                "max": 1,
                "scale": "Blues",
            },
            "restriction_population": {
                "title": "Restriction-QC population",
                "description": (
                    "Restriction-orientation fractions use post-dedup, "
                    "high-confidence UU pairs before contact filtering"
                ),
            },
            "overall_status": {
                "title": "Overall",
                "description": (
                    "PASS / FAIL / PASS_WITH_UNCERTAINTY / "
                    "PASS_WITH_NOT_ASSESSED"
                ),
            },
        },
        "data": {report["sample"]: fields},
    }


def differential_payload(report: dict) -> dict:
    """Expose differential inference eligibility rather than only its plots."""
    comparison = str(report["comparison"])
    fields = {
        "analysis_status": report["analysis_status"],
        "n_complete_pairs": int(report["n_complete_pairs"]),
        "publication_eligible": bool(report["publication_eligible"]),
        "publication_min_complete_pairs": int(
            report["publication_min_complete_pairs"]
        ),
        "design_formula": report.get("design_formula", ""),
        "candidate_loops": report.get("candidate_loops"),
        "tested_loops": report.get("tested_loops"),
    }
    return {
        "id": "differential_status",
        "section_name": "Differential analysis status",
        "description": (
            "Paired-design completeness and publication eligibility. "
            "PILOT_UNDERPOWERED results are exploratory even when model p-values "
            "are present."
        ),
        "plot_type": "table",
        "pconfig": {
            "id": "differential_status_table",
            "title": "Paired differential inference status",
        },
        "headers": {
            "analysis_status": {"title": "Analysis status"},
            "n_complete_pairs": {
                "title": "Complete pairs",
                "format": "{:,.0f}",
            },
            "publication_eligible": {"title": "Publication eligible"},
            "publication_min_complete_pairs": {
                "title": "Publication minimum",
                "format": "{:,.0f}",
            },
            "design_formula": {"title": "Fitted model"},
            "candidate_loops": {"title": "Candidate loops", "format": "{:,.0f}"},
            "tested_loops": {"title": "Tested loops", "format": "{:,.0f}"},
        },
        "data": {comparison: fields},
    }


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    report = json.loads(Path(snakemake.input.json).read_text())
    builders = {
        "apa": apa_payload,
        "balance": balance_payload,
        "differential": differential_payload,
        "loop_qc": loop_qc_payload,
    }
    try:
        builder = builders[snakemake.params.kind]
    except KeyError as exc:
        raise ValueError(f"Unknown MultiQC content kind: {snakemake.params.kind}") from exc
    payload = builder(report)

    output = Path(snakemake.output.json)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    log = Path(snakemake.log[0])
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(
        f"Wrote {snakemake.params.kind} custom content to {output}\n",
        encoding="utf-8",
    )


if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
