"""
Virtual 4C: pick a viewpoint (anchor of interest) and plot the contact frequency
of every other bin on the same chromosome from that viewpoint.
"""
from __future__ import annotations

import sys
from pathlib import Path

import cooler
import matplotlib.pyplot as plt
import numpy as np
import pyBigWig

sys.path.insert(0, str(Path(__file__).parent))
from utils import setup_logging  # noqa: E402


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])

    region = snakemake.params.region
    res = int(snakemake.params.res)
    chrom = region["chrom"]
    mid = (int(region["start"]) + int(region["end"])) // 2

    clr = cooler.Cooler(f"{snakemake.input.mcool}::resolutions/{res}")

    # Identify the viewpoint bin index
    bins = clr.bins().fetch(chrom).reset_index(drop=True)
    bin_row = bins[(bins["start"] <= mid) & (bins["end"] > mid)]
    if bin_row.empty:
        raise RuntimeError(f"No bin contains midpoint {mid} on {chrom}")
    bin_idx = int(bin_row.index[0])

    # Pull the row of contacts from the viewpoint to every bin on chrom
    mat = clr.matrix(balance=True).fetch(chrom)
    profile = np.nan_to_num(mat[bin_idx, :], nan=0.0).astype(float)
    starts = bins["start"].to_numpy(dtype=np.int64)
    ends = bins["end"].to_numpy(dtype=np.int64)

    # Write a bigWig of the v4C profile (chrom-restricted header is fine)
    bw = pyBigWig.open(str(snakemake.output.bw), "w")
    bw.addHeader([(chrom, int(clr.chromsizes[chrom]))])
    bw.addEntries(
        [chrom] * len(starts),
        starts.tolist(),
        ends=ends.tolist(),
        values=profile.tolist(),
    )
    bw.close()

    # Static figure restricted to the visualisation region
    fig, ax = plt.subplots(figsize=(10, 3))
    ax.fill_between(starts, profile, lw=0.5, color="#c0392b")
    ax.axvline(mid, color="k", ls="--", lw=0.5)
    ax.set_xlim(int(region["start"]), int(region["end"]))
    ax.set_xlabel(f"{chrom} position")
    ax.set_ylabel("normalised contacts")
    ax.set_title(f"Virtual 4C — {snakemake.wildcards.sample} @ {region['name']}")
    fig.tight_layout()
    fig.savefig(snakemake.output.png, dpi=150)


main(snakemake)  # type: ignore[name-defined]  # noqa: F821
