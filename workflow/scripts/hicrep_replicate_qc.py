"""
Stratum-adjusted Pearson correlation (HiCRep) between biological replicates.
HiCRep is the only matrix similarity metric robust to distance decay differences.

Inputs (from snakemake.input.mcools): one or more .mcool files for the same
subject + mark. Outputs (snakemake.output.json): per-pair SCC values.
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
from hicrep import hicrepSCC
from hicrep.utils import readMcool

sys.path.insert(0, str(Path(__file__).parent))
from utils import setup_logging, write_json  # noqa: E402


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    mcools = list(snakemake.input.mcools)
    bin_sz = int(snakemake.params.bin)
    max_dist = int(snakemake.params.maxd)
    h = int(snakemake.params.h)
    threshold = float(snakemake.config["hicrep"]["threshold_pass"])

    result = {
        "sample": snakemake.wildcards.sample,
        "bin_size": bin_sz,
        "max_dist": max_dist,
        "h": h,
        "n_replicates": len(mcools),
        "pairwise_scc": [],
        "mean_scc": None,
        "pass": False,
    }

    if len(mcools) < 2:
        # Single sample — no pairwise computation. Emit a placeholder.
        result["note"] = "Only one replicate; HiCRep skipped."
        write_json(result, snakemake.output.json)
        return

    sccs = []
    for a, b in combinations(mcools, 2):
        cool_a, _ = readMcool(a, bin_sz)
        cool_b, _ = readMcool(b, bin_sz)
        scc = hicrepSCC(cool_a, cool_b, h, max_dist, bin_sz)
        # hicrepSCC returns per-chromosome; mean across autosomes
        mean_scc = float(np.nanmean([v for k, v in scc.items() if "X" not in k and "Y" not in k]))
        result["pairwise_scc"].append({
            "a": Path(a).stem, "b": Path(b).stem, "scc": mean_scc
        })
        sccs.append(mean_scc)

    result["mean_scc"] = float(np.mean(sccs))
    result["pass"] = result["mean_scc"] >= threshold
    write_json(result, snakemake.output.json)


# Snakemake hands the `snakemake` object via global
main(snakemake)  # type: ignore[name-defined]  # noqa: F821
