"""
Stratum-adjusted Pearson correlation (HiCRep) between biological replicates.

Outputs a three-state QC result: PASS, FAIL, or NOT_ASSESSED. Single-replicate
samples are not treated as true passes.
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
        "threshold": threshold,
        "status": "NOT_ASSESSED",
        "pass": None,
    }

    if len(mcools) < 2:
        result["note"] = "Only one replicate; HiCRep was not assessed."
        write_json(result, snakemake.output.json)
        return

    sccs = []
    for a, b in combinations(mcools, 2):
        cool_a, _ = readMcool(a, bin_sz)
        cool_b, _ = readMcool(b, bin_sz)
        scc = hicrepSCC(cool_a, cool_b, h, max_dist, bin_sz)
        mean_scc = float(np.nanmean([v for k, v in scc.items() if "X" not in k and "Y" not in k]))
        result["pairwise_scc"].append({"a": Path(a).stem, "b": Path(b).stem, "scc": mean_scc})
        sccs.append(mean_scc)

    result["mean_scc"] = float(np.mean(sccs))
    result["pass"] = result["mean_scc"] >= threshold
    result["status"] = "PASS" if result["pass"] else "FAIL"
    write_json(result, snakemake.output.json)


main(snakemake)  # type: ignore[name-defined]  # noqa: F821
