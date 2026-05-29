"""
Export the E1 eigenvector (A/B compartments) from cooltools eigs-cis output
to a bigWig file. Pulled out of 06_loop_qc.smk to avoid fragile embedded
Python in shell.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pyBigWig

sys.path.insert(0, str(Path(__file__).parent))
from utils import read_chromsizes, setup_logging  # noqa: E402


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])

    df = pd.read_csv(snakemake.input.eigs, sep="\t").dropna(subset=["E1"])
    sizes = read_chromsizes(snakemake.params.chromsizes)
    # Drop rows whose chromosome is not in the chromsizes (Y / MT / decoys)
    df = df[df["chrom"].isin(sizes)]
    if df.empty:
        Path(snakemake.output.bw).write_bytes(b"")
        return
    bw = pyBigWig.open(str(snakemake.output.bw), "w")
    bw.addHeader([(c, int(s)) for c, s in sizes.items() if c in set(df["chrom"])])
    df = df.sort_values(["chrom", "start"])
    bw.addEntries(
        df["chrom"].astype(str).tolist(),
        df["start"].astype(int).tolist(),
        ends=df["end"].astype(int).tolist(),
        values=df["E1"].astype(float).tolist(),
    )
    bw.close()


main(snakemake)  # type: ignore[name-defined]  # noqa: F821
