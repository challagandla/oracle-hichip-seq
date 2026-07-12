"""
Aggregate Peak Analysis (APA) on a HiChIP loop set.

For each loop, extract a square of the contact matrix centred on the loop and
average across loops. The centre pixel is then compared against two backgrounds:

  apa_vs_random_shift  centre / centre of a random-shift control, where BOTH
                       anchors move by the same offset so the loop's genomic
                       separation is preserved. This is the score that decides
                       pass/fail -- it is distance-matched by construction.
  apa_score            centre / the two window corners that sit at the same
                       genomic separation as the centre. Reported for
                       comparability with published APA numbers.

Contact frequency inside the window is a function of distance from the diagonal:
a pixel (i, j) lies at separation D + (j - i) * bin_size. So corners at different
(j - i) are not interchangeable, and averaging all four -- as this did -- puts a
corner that is 2 * win * bin_size CLOSER to the diagonal into the denominator.
"""
from __future__ import annotations

import sys
from pathlib import Path

import cooler
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_loops_bedpe, setup_logging, write_json  # noqa: E402

MAX_SHIFT_BP = 1_000_000


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
        # The aggregate matrix is an output in its own right: the cohort figure
        # renders APA panels from it, and it must exist even when there is nothing
        # to aggregate, or one loopless sample takes the whole figure stage down.
        np.save(snakemake.output.npy, np.zeros((2 * win + 1, 2 * win + 1)))
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
            # Both anchors move by the SAME offset, which preserves the genomic
            # separation and so keeps the control on the same diagonal: the control
            # is distance-matched to the loop by construction. The offset must also
            # clear the window, or the "control" slides back onto the real loop and
            # the denominator absorbs the very signal it is meant to measure.
            shift = int(rng.integers(win_bp + bin_sz, MAX_SHIFT_BP))
            if rng.random() < 0.5:
                shift = -shift
            sq = _extract_square(clr, str(row.chrom1), mid1 + shift, mid2 + shift, win_bp)
            if sq is None or sq.shape != ctrl_agg.shape:
                continue
            ctrl_agg += np.nan_to_num(sq)
            n_ctrl_used += 1
        if n_ctrl_used:
            ctrl_agg /= n_ctrl_used
        ctrl_aggs.append(ctrl_agg)

    centre = agg[win, win]

    # Background must sit at the SAME genomic separation as the centre pixel.
    # In this window a pixel (i, j) lies at separation D + (j - i) * bin_sz, so the
    # four corners span D ± 2 * win * bin_sz -- at win=20 and 10 kb bins that is a
    # 400 kb swing. The bottom-left corner is therefore ~400 kb CLOSER to the
    # diagonal than the loop, carries far more contacts for purely distance-decay
    # reasons, and dominates a four-corner mean: the old denominator was inflated
    # and the APA score correspondingly crushed. Only the two corners on the
    # j - i = 0 anti-diagonal are distance-matched to the loop.
    background = np.concatenate([agg[:3, :3].ravel(), agg[-3:, -3:].ravel()])
    apa = float(centre / max(np.nanmean(background), 1e-9))

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

    # Pass/fail is decided on the random-shift control, not on the corner ratio.
    # The shifted control holds the loop's genomic separation fixed, so it is the
    # only one of the two that is distance-matched by construction; the corner
    # ratio is reported alongside it for comparability with published APA numbers.
    np.save(snakemake.output.npy, agg)

    score_min = float(snakemake.config["apa"]["score_min"])
    write_json({
        "sample": snakemake.wildcards.sample,
        "n_loops_used": int(n_used),
        "apa_score": apa,
        "apa_vs_random_shift": apa_vs_ctrl,
        "score_min": score_min,
        "pass": apa_vs_ctrl >= score_min,
    }, snakemake.output.json)


# Guarded so the module can be imported by the tests. Snakemake injects
# `snakemake` into the script's globals before executing it.
if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
