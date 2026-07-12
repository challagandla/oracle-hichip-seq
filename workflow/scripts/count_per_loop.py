"""
Count valid pairs supporting each union loop, per sample, from the unbalanced
cooler at the FitHiChIP bin size.

Two things this has to get right.

Anchors are not one bin wide. FitHiChIP runs with MergeInt=1, which merges
adjacent significant bins into one interaction, so an anchor is routinely 2-4
bins across. Counting only the bin containing `start` samples a fraction of each
loop's actual footprint, and that fraction is larger for wide anchors than for
narrow ones -- i.e. it tracks ChIP enrichment, which is exactly what differs
between the groups being compared, so it would enter DESeq2 as a
group-correlated bias. Each loop is therefore counted over the full rectangle of
bins its two anchors span.

The counting is done with one streaming pass over the pixel table rather than a
cooler .fetch() per loop. A union set here is O(10^5) loops and every .fetch() is
an indexed HDF5 range read; per-loop fetching costs hours per sample, whereas the
pixel table can be read linearly once.
"""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import cooler
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_loops_bedpe, setup_logging  # noqa: E402

CHUNK = 10_000_000  # pixel rows per read

log = logging.getLogger(__name__)


def _anchor_bins(clr: cooler.Cooler, chrom: str, start, end, res: int) -> np.ndarray:
    """Global bin ids spanned by [start, end) on `chrom`, clipped to the chromosome."""
    offset = int(clr.offset(chrom))
    n_chrom_bins = int(np.ceil(clr.chromsizes[chrom] / res))
    lo = max(int(start) // res, 0)
    hi = -(-int(end) // res)          # ceil division
    hi = max(hi, lo + 1)              # a zero-length anchor still covers its own bin
    hi = min(hi, n_chrom_bins)
    if hi <= lo:
        return np.empty(0, dtype=np.int64)
    return offset + np.arange(lo, hi, dtype=np.int64)


def count_loops(clr: cooler.Cooler, loops: pd.DataFrame, res: int) -> np.ndarray:
    """Pairs supporting each loop, summed over the full rectangle its anchors span."""
    if loops.empty:
        return np.zeros(0, dtype=np.int64)

    n_bins = int(clr.info["nbins"])
    known = set(clr.chromnames)

    # Expand every loop into the pixels its two anchors span. cooler stores the
    # upper triangle only, so each pair is emitted with bin1 <= bin2.
    pair_keys: list[np.ndarray] = []
    pair_loop: list[np.ndarray] = []
    for i, row in enumerate(loops.itertuples(index=False)):
        c1, c2 = str(row.chrom1), str(row.chrom2)
        if c1 not in known or c2 not in known:
            continue
        b1 = _anchor_bins(clr, c1, row.start1, row.end1, res)
        b2 = _anchor_bins(clr, c2, row.start2, row.end2, res)
        if b1.size == 0 or b2.size == 0:
            continue
        g1, g2 = np.meshgrid(b1, b2, indexing="ij")
        lo = np.minimum(g1, g2).ravel()
        hi = np.maximum(g1, g2).ravel()
        keys = np.unique(lo * n_bins + hi)  # unique: a self-overlapping anchor pair
        pair_keys.append(keys)              # must not count the same pixel twice
        pair_loop.append(np.full(keys.size, i, dtype=np.int64))

    counts = np.zeros(len(loops), dtype=np.int64)
    if not pair_keys:
        raise RuntimeError(
            "no loop anchor fell on the cooler bin grid -- check that the BEDPE and "
            "the cooler use the same chromosome naming and assembly"
        )

    # key -> loop is one-to-many: overlapping union anchors can share a pixel, and
    # each loop that claims it is credited with it.
    wanted = pd.DataFrame({"key": np.concatenate(pair_keys), "loop": np.concatenate(pair_loop)})
    wanted = wanted.sort_values("key", kind="stable").reset_index(drop=True)
    uniq_sorted = np.unique(wanted["key"].to_numpy())
    log.info("%d loops -> %d anchor pixels (%d distinct)",
             len(loops), len(wanted), uniq_sorted.size)

    pixels = clr.pixels()
    nnz = int(clr.info["nnz"])
    for lo_row in range(0, nnz, CHUNK):
        chunk = pixels[lo_row: lo_row + CHUNK]
        if chunk.empty:
            continue
        ck = (chunk["bin1_id"].to_numpy(dtype=np.int64) * n_bins
              + chunk["bin2_id"].to_numpy(dtype=np.int64))
        pos = np.minimum(np.searchsorted(uniq_sorted, ck), uniq_sorted.size - 1)
        hit = uniq_sorted[pos] == ck
        if not hit.any():
            continue
        matched = pd.DataFrame({"key": ck[hit], "count": chunk["count"].to_numpy()[hit]})
        agg = matched.merge(wanted, on="key", how="inner").groupby("loop")["count"].sum()
        counts[agg.index.to_numpy()] += agg.to_numpy(dtype=np.int64)

    return counts


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])

    res = int(snakemake.params.res)
    clr = cooler.Cooler(f"{snakemake.input.mcool}::resolutions/{res}")
    loops = load_loops_bedpe(snakemake.input.bedpe).reset_index(drop=True)

    out = loops.copy()
    if loops.empty:
        out["count"] = pd.Series(dtype="int64")
        out["sample"] = snakemake.wildcards.sample
        Path(snakemake.output.counts).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(snakemake.output.counts, sep="\t", index=False)
        log.warning("union BEDPE is empty; wrote an empty count table")
        return

    counts = count_loops(clr, loops, res)

    out["count"] = counts
    out["sample"] = snakemake.wildcards.sample
    log.info("sample=%s loops=%d nonzero=%d total_pairs=%d",
             snakemake.wildcards.sample, len(out),
             int((counts > 0).sum()), int(counts.sum()))
    Path(snakemake.output.counts).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(snakemake.output.counts, sep="\t", index=False)


# Guarded so the module can be imported by the tests. Snakemake injects `snakemake`
# into the script's globals before executing it; nothing else does.
if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
