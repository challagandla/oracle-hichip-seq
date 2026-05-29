"""
Aggregate Peak Analysis (APA) on a HiChIP loop set.

For each loop, extract a square of the contact matrix centred on the loop
and aggregate (mean) across all loops. Compare to N random-shift controls.
APA score = centre / mean(corners). Score ≥ 1.5 → good loops.
"""
from __future__ import annotations

import sys
from pathlib import Path

import cooler
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_loops_bedpe, setup_logging, write_json  # noqa: E402


def _extract_square(clr: cooler.Cooler, chrom: str, mid1: int, mid2: int,
                    win_bp: int, balanced: bool = True) -> np.ndarray | None:
    """Extract a (2W+1)×(2W+1) window of contact frequency around (mid1,mid2)."""
    try:
        region1 = (chrom, max(0, mid1 - win_bp), mid1 + win_bp)
        region2 = (chrom, max(0, mid2 - win_bp), mid2 + win_bp)
        mat = clr.matrix(balance=balanced).fetch(region1, region2)
    except Exception:
        return None
    return mat


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])

    bin_sz = int(snakemake.params.bin_size)
    win = int(snakemake.params.window)            # in bins
    win_bp = win * bin_sz
    min_dist = int(snakemake.params.min_dist)
    n_ctrl = int(snakemake.params.n_ctrl)

    clr = cooler.Cooler(f"{snakemake.input.mcool}::resolutions/{bin_sz}")
    loops = load_loops_bedpe(snakemake.input.loops)

    loops = loops[
        (loops.chrom1 == loops.chrom2) &
        ((loops.start2 - loops.start1).abs() >= min_dist)
    ].reset_index(drop=True)

    if len(loops) == 0:
        write_json({"sample": snakemake.wildcards.sample, "n_loops": 0,
                    "apa_score": None, "pass": False}, snakemake.output.json)
        plt.figure(); plt.title("No loops"); plt.savefig(snakemake.output.png); return

    # Aggregate real loops
    agg = np.zeros((2 * win + 1, 2 * win + 1))
    n_used = 0
    for _, row in loops.iterrows():
        mid1 = (int(row.start1) + int(row.end1)) // 2
        mid2 = (int(row.start2) + int(row.end2)) // 2
        sq = _extract_square(clr, str(row.chrom1), mid1, mid2, win_bp)
        if sq is None or sq.shape != (2 * win + 1, 2 * win + 1):
            continue
        sq = np.nan_to_num(sq, nan=0.0)
        agg += sq
        n_used += 1
    if n_used:
        agg /= n_used

    # Controls
    rng = np.random.default_rng(seed=42)
    ctrl_aggs = []
    for _ in range(n_ctrl):
        ctrl_agg = np.zeros_like(agg)
        n_ctrl_used = 0
        for _, row in loops.iterrows():
            mid1 = (int(row.start1) + int(row.end1)) // 2
            mid2 = (int(row.start2) + int(row.end2)) // 2
            shift = rng.integers(-1_000_000, 1_000_000)
            sq = _extract_square(clr, str(row.chrom1), mid1 + shift, mid2 + shift, win_bp)
            if sq is None or sq.shape != ctrl_agg.shape:
                continue
            ctrl_agg += np.nan_to_num(sq)
            n_ctrl_used += 1
        if n_ctrl_used:
            ctrl_agg /= n_ctrl_used
        ctrl_aggs.append(ctrl_agg)

    centre = agg[win, win]
    corners = np.concatenate([
        agg[:3, :3].ravel(), agg[:3, -3:].ravel(),
        agg[-3:, :3].ravel(), agg[-3:, -3:].ravel(),
    ])
    apa = float(centre / max(np.nanmean(corners), 1e-9))

    ctrl_centres = [c[win, win] for c in ctrl_aggs if c.any()]
    apa_vs_ctrl = float(centre / max(np.nanmean(ctrl_centres), 1e-9))

    # Plot
    fig, ax = plt.subplots(figsize=(4, 4))
    im = ax.imshow(np.log2(agg + 1), cmap="Reds", origin="lower",
                   extent=[-win, win, -win, win])
    ax.set_title(f"APA {snakemake.wildcards.sample}\nscore={apa:.2f}, vs_ctrl={apa_vs_ctrl:.2f}, n={n_used}")
    ax.set_xlabel(f"bins ({bin_sz//1000} kb)"); ax.set_ylabel("bins")
    fig.colorbar(im, ax=ax, label="log2(1+contacts)")
    fig.tight_layout()
    fig.savefig(snakemake.output.png, dpi=150)

    write_json({
        "sample": snakemake.wildcards.sample,
        "n_loops_used": int(n_used),
        "apa_score": apa,
        "apa_vs_random_shift": apa_vs_ctrl,
        "pass": apa >= 1.5,
    }, snakemake.output.json)


main(snakemake)  # type: ignore[name-defined]  # noqa: F821
