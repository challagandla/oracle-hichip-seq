"""
Count valid pairs supporting each union loop, per sample, from the unbalanced
cooler at the FitHiChIP bin size.

Each loop is counted over the full rectangle of bins its two anchors span, not over
the single bin containing `start`.

FitHiChIP 11.0 normally emits one-bin anchors, including in its merged output.

The rectangle is kept anyway, because the union BEDPE is not required to come from
FitHiChIP. Any caller that emits a wider anchor -- mustache, a lifted-over published
loop set, a merged consensus across resolutions -- would otherwise be counted on an
arbitrary fraction of its own footprint, and that fraction scales with anchor width,
which scales with ChIP enrichment, which is what differs between the groups being
compared. That is a group-correlated bias rather than noise, and it is silent.

Counting uses one streaming pass over the pixel table instead of a cooler .fetch()
per loop. A union set is O(10^5) loops and each
.fetch() is an indexed HDF5 range read, so per-loop fetching costs hours per sample;
the pixel table is read linearly once.
"""
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


def count_loops(
    clr: cooler.Cooler,
    loops: pd.DataFrame,
    res: int,
    lower_distance: int,
    upper_distance: int,
) -> np.ndarray:
    """Count supporting pairs inside each footprint and the caller distance range.

    Tolerance-expanded anchors can straddle FitHiChIP's lower or upper distance
    boundary.  Pixels outside that exact search range were never eligible for the
    caller's background model and must not enter the differential count matrix.
    """
    if not (0 <= int(lower_distance) < int(upper_distance)):
        raise ValueError("distance bounds must satisfy 0 <= lower < upper")
    if loops.empty:
        return np.zeros(0, dtype=np.int64)

    n_bins = int(clr.info["nbins"])
    known = set(clr.chromnames)
    bin_starts = clr.bins()[:]["start"].to_numpy(dtype=np.int64)

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
        distance = np.abs(bin_starts[hi] - bin_starts[lo])
        eligible = (
            (distance >= int(lower_distance))
            & (distance <= int(upper_distance))
        )
        lo = lo[eligible]
        hi = hi[eligible]
        if lo.size == 0:
            continue
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

    counts = count_loops(
        clr,
        loops,
        res,
        int(snakemake.params.lower_distance),
        int(snakemake.params.upper_distance),
    )

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
