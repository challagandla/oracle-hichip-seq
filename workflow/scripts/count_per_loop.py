"""
For a union BEDPE of loops and a sample's .mcool, count valid pairs supporting
each loop anchor pair at the configured resolution.
"""
from __future__ import annotations

import sys
from pathlib import Path

import cooler
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_loops_bedpe, setup_logging  # noqa: E402


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])

    res = int(snakemake.params.res)
    clr = cooler.Cooler(f"{snakemake.input.mcool}::resolutions/{res}")
    loops = load_loops_bedpe(snakemake.input.bedpe)

    counts: list[int] = []
    for _, row in loops.iterrows():
        try:
            mat = clr.matrix(balance=False).fetch(
                (str(row.chrom1), int(row.start1), int(row.end1)),
                (str(row.chrom2), int(row.start2), int(row.end2)),
            )
            counts.append(int(mat.sum()))
        except Exception:
            counts.append(0)

    out = loops.copy()
    out["count"] = counts
    out["sample"] = snakemake.wildcards.sample
    Path(snakemake.output.counts).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(snakemake.output.counts, sep="\t", index=False)


main(snakemake)  # type: ignore[name-defined]  # noqa: F821
