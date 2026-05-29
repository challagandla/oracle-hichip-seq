"""
Differential loop analysis between cases (e.g. tumor) and controls (e.g. healthy).

Default: pyDESeq2 on the loop-by-sample count matrix.
Alternative: invoke diffHic via Rscript (requires bioconductor-diffhic installed).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils import setup_logging  # noqa: E402


def _load_count_table(files: list[str]) -> tuple[pd.DataFrame, list[str]]:
    """Read all per-sample counts.tsv into a single (loop × sample) matrix."""
    frames = []
    sample_order: list[str] = []
    for f in files:
        df = pd.read_csv(f, sep="\t")
        sid = df["sample"].iloc[0]
        sample_order.append(sid)
        key = df[["chrom1", "start1", "end1", "chrom2", "start2", "end2"]].astype(str).agg("_".join, axis=1)
        frames.append(pd.Series(df["count"].values, index=key, name=sid))
    M = pd.concat(frames, axis=1).fillna(0).astype(int)
    return M, sample_order


def _pydeseq2(M: pd.DataFrame, cases: list[str], controls: list[str], fdr: float, lfc: float) -> pd.DataFrame:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    metadata = pd.DataFrame(
        {"condition": ["case" if s in cases else "control" for s in M.columns]},
        index=M.columns,
    )
    counts = M.T  # samples in rows for pyDESeq2
    dds = DeseqDataSet(counts=counts, metadata=metadata, design_factors="condition", quiet=True)
    dds.deseq2()
    stat = DeseqStats(dds, contrast=("condition", "case", "control"), quiet=True)
    stat.summary()
    res = stat.results_df.copy()
    res["sig"] = (res["padj"] < fdr) & (res["log2FoldChange"].abs() >= lfc)
    return res


def _volcano(res: pd.DataFrame, out_png: str | Path, fdr: float, lfc: float) -> None:
    fig, ax = plt.subplots(figsize=(5, 4.5))
    x = res["log2FoldChange"]
    y = -np.log10(res["padj"].replace(0, np.nextafter(0, 1)))
    sig = res["sig"]
    ax.scatter(x[~sig], y[~sig], s=4, c="#ccc", alpha=0.5)
    ax.scatter(x[sig & (x > 0)], y[sig & (x > 0)], s=6, c="#c0392b", label=f"up (n={int((sig & (x>0)).sum())})")
    ax.scatter(x[sig & (x < 0)], y[sig & (x < 0)], s=6, c="#1f5fbf", label=f"down (n={int((sig & (x<0)).sum())})")
    ax.axhline(-np.log10(fdr), c="k", ls="--", lw=0.5)
    ax.axvline( lfc, c="k", ls="--", lw=0.5); ax.axvline(-lfc, c="k", ls="--", lw=0.5)
    ax.set_xlabel("log2 fold change (case / control)")
    ax.set_ylabel("-log10 adjusted p")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)


def _ma(res: pd.DataFrame, out_png: str | Path) -> None:
    fig, ax = plt.subplots(figsize=(5, 4))
    base = np.log10(res["baseMean"].replace(0, np.nextafter(0, 1)))
    ax.scatter(base, res["log2FoldChange"], s=4, c=np.where(res["sig"], "#c0392b", "#bbb"), alpha=0.6)
    ax.axhline(0, c="k", lw=0.5)
    ax.set_xlabel("log10(base mean)"); ax.set_ylabel("log2 FC")
    fig.tight_layout(); fig.savefig(out_png, dpi=150)


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])

    method = snakemake.params.method
    fdr = float(snakemake.params.fdr)
    lfc = float(snakemake.params.log2fc_min)
    cases = list(snakemake.params.cases)
    controls = list(snakemake.params.controls)

    M, order = _load_count_table(list(snakemake.input.counts))

    if method == "pyDESeq2":
        res = _pydeseq2(M, cases, controls, fdr=fdr, lfc=lfc)
    elif method == "diffhic_r":
        raise NotImplementedError("diffhic_r path: call Rscript scripts/differential_diffhic.R")
    else:
        raise ValueError(f"Unknown differential method {method!r}")

    # Decorate result with loop coordinates
    loop_keys = M.index.str.split("_", expand=True)
    res = res.reset_index().rename(columns={"index": "loop_key"})
    coords = pd.DataFrame(loop_keys.tolist(),
                          columns=["chrom1", "start1", "end1", "chrom2", "start2", "end2"],
                          index=M.index)
    res = res.join(coords, on="loop_key").sort_values("padj")
    res.to_csv(snakemake.output.tsv, sep="\t", index=False)

    _volcano(res, snakemake.output.volcano, fdr=fdr, lfc=lfc)
    _ma(res, snakemake.output.ma)


main(snakemake)  # type: ignore[name-defined]  # noqa: F821
