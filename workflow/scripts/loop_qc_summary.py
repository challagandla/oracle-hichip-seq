"""
Aggregate per-sample QC into a single JSON + Markdown summary with pass/fail
flags. Replicate QC uses PASS / FAIL / NOT_ASSESSED instead of pretending a
single-replicate sample passed HiCRep.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils import setup_logging, write_json  # noqa: E402


TH = {
    "valid_pair_yield_pct": 25.0,
    "duplicate_pct_max": 50.0,
    "cis_fraction_min": 0.70,
    "apa_score_min": 1.5,
    # hicrep_scc_min is NOT hardcoded here: it is taken from config at runtime.
    # A second copy of the threshold in this file silently overrode the configured
    # one in the report, so the two could disagree about whether a sample passed.
    "n_loops_min": 1000,
}


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


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])

    TH["hicrep_scc_min"] = float(snakemake.config["hicrep"]["threshold_pass"])

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
    apa_score = apa.get("apa_score")
    hicrep_mean = hicrep.get("mean_scc")

    pass_flags = {
        "valid_pair_yield": valid_yield >= TH["valid_pair_yield_pct"],
        "duplicate_pct": dup_pct <= TH["duplicate_pct_max"],
        "cis_fraction": cis_frac >= TH["cis_fraction_min"],
        "n_loops": n_loops >= TH["n_loops_min"],
        "apa_score": (apa_score or 0) >= TH["apa_score_min"],
    }
    status_flags = {k: _status(v) for k, v in pass_flags.items()}
    if hicrep_mean is None:
        status_flags["hicrep_scc"] = "NOT_ASSESSED"
        hicrep_pass = None
    else:
        hicrep_pass = hicrep_mean >= TH["hicrep_scc_min"]
        status_flags["hicrep_scc"] = _status(hicrep_pass)

    hard_pass = all(pass_flags.values()) and (hicrep_pass is not False)
    if not hard_pass:
        overall_status = "FAIL"
    elif "NOT_ASSESSED" in status_flags.values():
        overall_status = "PASS_WITH_NOT_ASSESSED"
    else:
        overall_status = "PASS"

    report = {
        "sample": snakemake.wildcards.sample,
        "valid_pair_yield_pct": float(valid_yield),
        "duplicate_pct": float(dup_pct),
        "cis_fraction": float(cis_frac),
        "n_loops": int(n_loops),
        "apa_score": apa_score,
        "apa_vs_random_shift": apa.get("apa_vs_random_shift"),
        "hicrep_mean_scc": hicrep_mean,
        "hicrep_status": hicrep.get("status", status_flags["hicrep_scc"]),
        "thresholds": TH,
        "pass_flags": pass_flags,
        "status_flags": status_flags,
        "overall_status": overall_status,
        "overall_pass": overall_status.startswith("PASS"),
    }
    write_json(report, snakemake.output.json)

    apa_line = f"- APA score: **{apa_score:.2f}** (threshold ≥ {TH['apa_score_min']})" if apa_score is not None else "- APA score: **NA**"
    hicrep_line = f"- HiCRep mean SCC: **{hicrep_mean:.3f}** (threshold ≥ {TH['hicrep_scc_min']})" if hicrep_mean is not None else "- HiCRep mean SCC: **NOT_ASSESSED**"
    md_lines = [
        f"# QC report — {snakemake.wildcards.sample}",
        "",
        f"- Valid pair yield: **{valid_yield:.1f}%** (threshold ≥ {TH['valid_pair_yield_pct']}%)",
        f"- Duplicate %: **{dup_pct:.1f}%** (threshold ≤ {TH['duplicate_pct_max']}%)",
        f"- Cis fraction: **{cis_frac:.2f}** (threshold ≥ {TH['cis_fraction_min']:.2f})",
        f"- N significant loops: **{n_loops}** (threshold ≥ {TH['n_loops_min']})",
        apa_line,
        hicrep_line,
        "",
        f"**Overall: {overall_status}**",
    ]
    Path(snakemake.output.md).write_text("\n".join(md_lines))


main(snakemake)  # type: ignore[name-defined]  # noqa: F821
