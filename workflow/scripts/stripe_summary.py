#!/usr/bin/env python3
"""Aggregate per-sample stripe calls into one table.

Reported per sample, alongside the anchor mark, because stripe counts are only
comparable within an anchor type: an architectural anchor (CTCF, cohesin) holds
the extruding complex and produces stripes by definition, whereas an H3K27ac
anchor set marks enhancers, which are not extrusion anchors and yield far fewer.
A lower stripe count on H3K27ac is the expected result, not a worse experiment.
"""
from __future__ import annotations

import pandas as pd

samples = list(snakemake.params.samples)  # noqa: F821
marks = dict(snakemake.params.marks)  # noqa: F821
res = int(snakemake.params.res)  # noqa: F821

rows = []
for path, sample in zip(snakemake.input, samples):  # noqa: F821
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        df = pd.DataFrame()

    n = len(df)
    length_col = "length" if "length" in df.columns else None
    rows.append(
        {
            "sample": sample,
            "mark": marks.get(sample, "NA"),
            "resolution": res,
            "n_stripes": n,
            "median_length_kb": (
                round(float(df[length_col].median()) / 1000, 1)
                if n and length_col
                else 0.0
            ),
            "max_length_kb": (
                round(float(df[length_col].max()) / 1000, 1)
                if n and length_col
                else 0.0
            ),
        }
    )

out = pd.DataFrame(rows).sort_values(["mark", "sample"])
out.to_csv(snakemake.output.tsv, sep="\t", index=False)  # noqa: F821

with open(snakemake.log[0], "w") as fh:  # noqa: F821
    fh.write(out.to_string(index=False) + "\n")
    arch = out[out["mark"].isin(["CTCF", "SMC1A", "RAD21"])]
    hist = out[~out["mark"].isin(["CTCF", "SMC1A", "RAD21"])]
    if len(arch) and len(hist):
        fh.write(
            f"\narchitectural anchors: median {arch.n_stripes.median():.0f} stripes; "
            f"histone anchors: median {hist.n_stripes.median():.0f}. "
            "Fewer stripes on a histone anchor set is expected.\n"
        )
