"""Export the E1 eigenvector (A/B compartments) to bigWig."""
import sys
from pathlib import Path

import pandas as pd
import pyBigWig

sys.path.insert(0, str(Path(__file__).parent))
from utils import read_chromsizes, setup_logging  # noqa: E402


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    out = Path(snakemake.output.bw)
    out.parent.mkdir(parents=True, exist_ok=True)

    sizes = read_chromsizes(snakemake.params.chromsizes)
    df = pd.read_csv(snakemake.input.eigs, sep="\t")
    if "E1" not in df.columns:
        raise ValueError(f"{snakemake.input.eigs} does not contain an E1 column")
    df = df.dropna(subset=["E1"])
    df = df[df["chrom"].isin(sizes)]

    bw = pyBigWig.open(str(out), "w")
    header_chroms = sorted(set(df["chrom"])) if not df.empty else list(sizes.keys())
    bw.addHeader([(c, int(sizes[c])) for c in header_chroms if c in sizes])
    if not df.empty:
        df = df.sort_values(["chrom", "start"])
        bw.addEntries(
            df["chrom"].astype(str).tolist(),
            df["start"].astype(int).tolist(),
            ends=df["end"].astype(int).tolist(),
            values=df["E1"].astype(float).tolist(),
        )
    bw.close()


# Guarded so the module can be imported by the tests. Snakemake injects
# `snakemake` into the script's globals before executing it.
if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
