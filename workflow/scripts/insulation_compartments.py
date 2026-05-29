"""
Standalone helper to plot insulation profile + A/B compartment track for a
sample across a genome region. Useful for QC review notebooks.
Not invoked by the main Snakefile (those are produced by cooltools rules).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils import setup_logging  # noqa: E402


def plot_region(insul_tsv: str | Path, eigs_tsv: str | Path,
                chrom: str, start: int, end: int, out_png: str | Path) -> None:
    insul = pd.read_csv(insul_tsv, sep="\t")
    eigs = pd.read_csv(eigs_tsv, sep="\t")

    sub_i = insul[(insul.chrom == chrom) & (insul.start >= start) & (insul.end <= end)]
    sub_e = eigs[(eigs.chrom == chrom) & (eigs.start >= start) & (eigs.end <= end)]

    fig, axes = plt.subplots(2, 1, figsize=(10, 4), sharex=True)
    ins_col = "log2_insulation_score" if "log2_insulation_score" in sub_i.columns else sub_i.columns[-1]
    axes[0].plot(sub_i.start, sub_i[ins_col], lw=0.8, color="#444")
    axes[0].set_ylabel(ins_col); axes[0].axhline(0, c="grey", lw=0.5)

    axes[1].fill_between(sub_e.start, 0, sub_e.E1,
                         where=(sub_e.E1 > 0), color="#c0392b", alpha=0.6)
    axes[1].fill_between(sub_e.start, 0, sub_e.E1,
                         where=(sub_e.E1 <= 0), color="#1f5fbf", alpha=0.6)
    axes[1].set_ylabel("E1 (A/B)"); axes[1].set_xlabel(f"{chrom} position")

    fig.tight_layout()
    fig.savefig(out_png, dpi=150)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--insul", required=True)
    p.add_argument("--eigs", required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("--start", type=int, required=True)
    p.add_argument("--end", type=int, required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()
    setup_logging(None)
    plot_region(args.insul, args.eigs, args.chrom, args.start, args.end, args.out)
