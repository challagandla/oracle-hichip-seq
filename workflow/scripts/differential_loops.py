"""
Differential loop analysis between explicitly configured case/control groups.

Default: pyDESeq2 on the loop-by-sample count matrix. The Snakemake rule guards
comparison definition so groups are mark/tissue/protocol compatible before this
script is called.
"""
from __future__ import annotations

import json
import hashlib
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils import setup_logging  # noqa: E402

LOOP_COORD_COLS = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]


def _loop_keys(coords: pd.DataFrame) -> pd.Series:
    encoded = coords[LOOP_COORD_COLS].astype(str).agg("\x1f".join, axis=1)
    return encoded.map(lambda value: hashlib.sha1(value.encode("utf-8")).hexdigest())


def _load_count_table(files: list[str]) -> tuple[pd.DataFrame, list[str], pd.DataFrame]:
    """Read all per-sample counts.tsv into a single (loop × sample) matrix."""
    frames = []
    coord_frames = []
    sample_order: list[str] = []
    for f in files:
        df = pd.read_csv(f, sep="\t")
        if df.empty:
            continue
        sid = df["sample"].iloc[0]
        sample_order.append(sid)
        coords = df[LOOP_COORD_COLS].copy()
        key = _loop_keys(coords)
        frames.append(pd.Series(df["count"].values, index=key, name=sid))
        coords.insert(0, "loop_key", key)
        coord_frames.append(coords)
    if not frames:
        raise RuntimeError("No loop-count tables contained data")
    M = pd.concat(frames, axis=1).fillna(0).astype(int)
    coords = pd.concat(coord_frames, ignore_index=True).drop_duplicates("loop_key")
    return M, sample_order, coords


def _pydeseq2(M: pd.DataFrame, cases: list[str], controls: list[str], fdr: float, lfc: float) -> pd.DataFrame:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    metadata = pd.DataFrame(
        {"condition": ["case" if s in cases else "control" for s in M.columns]},
        index=M.columns,
    )
    if metadata["condition"].nunique() != 2:
        raise RuntimeError("Differential loop analysis requires both case and control samples in the count matrix")
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
    ax.axvline(lfc, c="k", ls="--", lw=0.5)
    ax.axvline(-lfc, c="k", ls="--", lw=0.5)
    ax.set_xlabel("log2 fold change (case / control)")
    ax.set_ylabel("-log10 adjusted p")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)


def _ma(res: pd.DataFrame, out_png: str | Path) -> None:
    fig, ax = plt.subplots(figsize=(5, 4))
    base = np.log10(res["baseMean"].replace(0, np.nextafter(0, 1)))
    ax.scatter(base, res["log2FoldChange"], s=4, c=np.where(res["sig"], "#c0392b", "#bbb"), alpha=0.6)
    ax.axhline(0, c="k", lw=0.5)
    ax.set_xlabel("log10(base mean)")
    ax.set_ylabel("log2 FC")
    fig.tight_layout()
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])

    method = snakemake.params.method
    fdr = float(snakemake.params.fdr)
    lfc = float(snakemake.params.log2fc_min)
    cases, controls, comp = snakemake.params.groups
    cases = list(cases)
    controls = list(controls)
    comp = dict(comp)

    if not cases or not controls:
        raise RuntimeError(f"Comparison {snakemake.wildcards.comparison} has empty case/control groups")

    M, order, coords = _load_count_table(list(snakemake.input.counts))
    missing_cases = sorted(set(cases) - set(M.columns))
    missing_controls = sorted(set(controls) - set(M.columns))
    if missing_cases or missing_controls:
        raise RuntimeError(f"Missing loop-count columns: cases={missing_cases}; controls={missing_controls}")

    if method == "pyDESeq2":
        res = _pydeseq2(M, cases, controls, fdr=fdr, lfc=lfc)
    elif method == "diffhic_r":
        raise NotImplementedError("diffhic_r path: call Rscript scripts/differential_diffhic.R")
    else:
        raise ValueError(f"Unknown differential method {method!r}")

    res.index = res.index.astype(str)
    res.index.name = "loop_key"
    res = res.reset_index().merge(coords, on="loop_key", how="left").sort_values("padj")
    Path(snakemake.output.tsv).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(snakemake.output.tsv, sep="\t", index=False)

    design = {"comparison": snakemake.wildcards.comparison, "cases": cases, "controls": controls, "config": comp}
    Path(snakemake.output.design).write_text(json.dumps(design, indent=2))
    _volcano(res, snakemake.output.volcano, fdr=fdr, lfc=lfc)
    _ma(res, snakemake.output.ma)


main(snakemake)  # type: ignore[name-defined]  # noqa: F821
