"""
Virtual 4C: pick a viewpoint (anchor of interest) and plot the contact frequency
of every other bin on the same chromosome from that viewpoint.
"""
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

    # Identify the viewpoint bin
    bins = clr.bins().fetch(chrom).reset_index(drop=True)
    bin_row = bins[(bins["start"] <= mid) & (bins["end"] > mid)]
    if bin_row.empty:
        raise RuntimeError(f"No bin contains midpoint {mid} on {chrom}")
    vp_start = int(bin_row["start"].iloc[0])
    vp_end = int(bin_row["end"].iloc[0])

    # Fetch ONLY the viewpoint row against the chromosome, never the whole
    # chromosome. fetch(chrom) returns a dense n x n array: chr1 at 5 kb is ~49,700
    # bins, so that single call asks for ~20 GB of RAM to read one row out of it,
    # and several viz jobs run concurrently.
    row = clr.matrix(balance=True).fetch((chrom, vp_start, vp_end), chrom)
    profile = np.nan_to_num(np.asarray(row, dtype=float).ravel(), nan=0.0)
    starts = bins["start"].to_numpy(dtype=np.int64)
    ends = bins["end"].to_numpy(dtype=np.int64)
    if profile.size != starts.size:
        raise RuntimeError(
            f"virtual 4C profile has {profile.size} values for {starts.size} bins on {chrom}"
        )

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


# Guarded so the module can be imported by the tests. Snakemake injects
# `snakemake` into the script's globals before executing it.
if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
