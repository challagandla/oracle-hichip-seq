"""
Aggregate per-sample QC into a single JSON + Markdown summary with pass/fail flags.
Cutoffs follow the BEST_PRACTICES doc.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils import setup_logging, write_json  # noqa: E402


# ----- thresholds (mirror docs/BEST_PRACTICES.md) -----
TH = {
    "valid_pair_yield_pct": 25.0,     # ≥ 25% of input read pairs end up as UU
    "duplicate_pct_max": 50.0,        # > 50% dedup → undersampled library
    "cis_fraction_min": 0.70,         # ≥ 70% cis reads
    "apa_score_min": 1.5,
    "hicrep_scc_min": 0.85,
    "n_loops_min": 1000,              # ≥ 1k FDR-significant loops per HiChIP sample
}


def parse_pairtools_stats(path: str | Path) -> dict[str, float]:
    """Parse the key=value style file emitted by `pairtools stats`."""
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


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])

    pair_stats = parse_pairtools_stats(snakemake.input.pair_stats)
    dedup_stats = parse_pairtools_stats(snakemake.input.dedup_stats)
    apa = json.loads(Path(snakemake.input.apa_json).read_text())
    hicrep = json.loads(Path(snakemake.input.hicrep).read_text())
    loops = pd.read_csv(snakemake.input.loops, sep="\t", header=None, comment="#")

    total = pair_stats.get("total", 0.0)
    total_dups = dedup_stats.get("total_dups", 0.0)
    cis = pair_stats.get("cis", 0.0)
    n_loops = len(loops)

    valid_yield = (pair_stats.get("total_nodups", 0.0) / total * 100.0) if total else 0.0
    dup_pct = (total_dups / max(total, 1)) * 100.0
    cis_frac = cis / max(pair_stats.get("total_nodups", 1), 1)

    report = {
        "sample": snakemake.wildcards.sample,
        "valid_pair_yield_pct": float(valid_yield),
        "duplicate_pct": float(dup_pct),
        "cis_fraction": float(cis_frac),
        "n_loops": int(n_loops),
        "apa_score": apa.get("apa_score"),
        "apa_vs_random_shift": apa.get("apa_vs_random_shift"),
        "hicrep_mean_scc": hicrep.get("mean_scc"),
        "thresholds": TH,
        "pass_flags": {
            "valid_pair_yield": valid_yield >= TH["valid_pair_yield_pct"],
            "duplicate_pct":    dup_pct    <= TH["duplicate_pct_max"],
            "cis_fraction":     cis_frac   >= TH["cis_fraction_min"],
            "n_loops":          n_loops    >= TH["n_loops_min"],
            "apa_score":        (apa.get("apa_score") or 0)   >= TH["apa_score_min"],
            "hicrep_scc":       (hicrep.get("mean_scc") or 1) >= TH["hicrep_scc_min"]
                                  if hicrep.get("mean_scc") is not None else True,
        },
    }
    report["overall_pass"] = all(report["pass_flags"].values())

    write_json(report, snakemake.output.json)

    md_lines = [
        f"# QC report — {snakemake.wildcards.sample}",
        "",
        f"- Valid pair yield: **{valid_yield:.1f}%** (threshold ≥ {TH['valid_pair_yield_pct']}%)",
        f"- Duplicate %: **{dup_pct:.1f}%** (threshold ≤ {TH['duplicate_pct_max']}%)",
        f"- Cis fraction: **{cis_frac:.2f}** (threshold ≥ {TH['cis_fraction_min']:.2f})",
        f"- N significant loops: **{n_loops}** (threshold ≥ {TH['n_loops_min']})",
        f"- APA score: **{apa.get('apa_score'):.2f}** (threshold ≥ {TH['apa_score_min']})"
            if apa.get('apa_score') is not None else "- APA score: NA",
        f"- HiCRep mean SCC: **{hicrep.get('mean_scc')}** (threshold ≥ {TH['hicrep_scc_min']})",
        "",
        f"**Overall: {'PASS' if report['overall_pass'] else 'FAIL'}**",
    ]
    Path(snakemake.output.md).write_text("\n".join(md_lines))


main(snakemake)  # type: ignore[name-defined]  # noqa: F821
