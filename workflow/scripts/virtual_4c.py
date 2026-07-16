"""Plot a chromosome-wide virtual-4C profile from an explicit viewpoint.

The viewpoint is configured separately from the plotting window. A region's
midpoint is rarely a biologically meaningful anchor and must never be substituted
for its promoter or enhancer coordinate implicitly.
"""
import sys
from pathlib import Path

import cooler
import matplotlib.pyplot as plt
import numpy as np
import pyBigWig

sys.path.insert(0, str(Path(__file__).parent))
from balance_utils import load_balance_report, resolution_balance  # noqa: E402
from utils import setup_logging  # noqa: E402
from viz_viewpoint import resolve_viewpoint  # noqa: E402


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])

    region = snakemake.params.region
    res = int(snakemake.params.res)
    chrom = region["chrom"]
    viewpoint, viewpoint_label = resolve_viewpoint(region)
    balance = resolution_balance(
        load_balance_report(snakemake.input.balance), res
    )

    clr = cooler.Cooler(f"{snakemake.input.mcool}::resolutions/{res}")
    if chrom not in clr.chromnames:
        raise ValueError(f"chromosome {chrom!r} is absent from {snakemake.input.mcool}")
    if viewpoint >= int(clr.chromsizes[chrom]):
        raise ValueError(
            f"viewpoint {chrom}:{viewpoint} exceeds chromosome length "
            f"{int(clr.chromsizes[chrom])}"
        )

    # Identify the viewpoint bin
    bins = clr.bins().fetch(chrom).reset_index(drop=True)
    bin_row = bins[(bins["start"] <= viewpoint) & (bins["end"] > viewpoint)]
    if bin_row.empty:
        raise RuntimeError(f"no {res:,}-bp bin contains viewpoint {chrom}:{viewpoint}")
    vp_start = int(bin_row["start"].iloc[0])
    vp_end = int(bin_row["end"].iloc[0])

    # Fetch ONLY the viewpoint row against the chromosome, never the whole
    # chromosome. fetch(chrom) returns a dense n x n array: chr1 at 5 kb is ~49,700
    # bins, so that single call asks for ~20 GB of RAM to read one row out of it,
    # and several viz jobs run concurrently.
    row = clr.matrix(balance=bool(balance["use_balanced"])).fetch(
        (chrom, vp_start, vp_end), chrom
    )
    profile = np.nan_to_num(
        np.asarray(row, dtype=float).ravel(), nan=0.0, posinf=0.0, neginf=0.0
    )
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
    ax.axvline(viewpoint, color="k", ls="--", lw=0.7, label=viewpoint_label)
    ax.set_xlim(int(region["start"]), int(region["end"]))
    ax.set_xlabel(f"{chrom} position")
    ax.set_ylabel(
        "ICE-balanced contact signal"
        if balance["use_balanced"] else "raw contact count"
    )
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    ax.set_title(
        f"Virtual 4C — {snakemake.wildcards.sample} from "
        f"{viewpoint_label} ({chrom}:{viewpoint:,}; {res // 1000} kb bin)\n"
        f"{balance['normalization']} (balance {balance['status']})"
    )
    fig.tight_layout()
    fig.savefig(snakemake.output.png, dpi=150)


# Guarded so the module can be imported by the tests. Snakemake injects
# `snakemake` into the script's globals before executing it.
if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
