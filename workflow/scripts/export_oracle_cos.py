"""
Export a HiChIP sample to the versioned ORACLE Chromatin Operating System
(COS) HiChIP prototype representation.

Current scope: HiChIP structural graph + peak-overlap prototype node features.
True continuous per-mark signal tracks (bigWig/RPKM), ATAC and methylation
channels are intentionally left for compatible modality-specific exports and later merge.
"""
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_loops_bedpe, read_chromsizes, select_insulation_column, setup_logging  # noqa: E402


NODE_FEATURE_CHANNELS = [
    "peak_overlap_count_per_kb", "insulation", "E1_eigenvector",
]
PRIMARY_CHROMOSOME = re.compile(r"^chr(?:[1-9]|1[0-9]|2[0-2]|X)$")


def _sha256(path: str | Path) -> str:
    """Return hex SHA-256 of a file, or ``MISSING`` when absent."""
    p = Path(path)
    if not p.exists():
        return "MISSING"
    digest = hashlib.sha256()
    with p.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _single_annotation(table: pd.DataFrame, column: str, default=None):
    if column not in table.columns:
        return default
    values = table[column].dropna().astype(str).unique().tolist()
    if not values:
        return default
    if len(values) != 1:
        raise ValueError(f"{column} has inconsistent values: {values}")
    return values[0]


def node_feature_availability(
    insulation: pd.DataFrame,
    insulation_value_column: str,
    eigs: pd.DataFrame,
    eigs_status: dict,
) -> dict:
    """Describe whether each unchanged numeric channel is interpretable."""
    insulation_values = pd.to_numeric(
        insulation.get(insulation_value_column, pd.Series(dtype=float)),
        errors="coerce",
    )
    insulation_available = bool(np.isfinite(insulation_values).any())
    e1_values = pd.to_numeric(eigs.get("E1", pd.Series(dtype=float)), errors="coerce")
    e1_available = bool(eigs_status.get("available")) and bool(
        np.isfinite(e1_values).any()
    )
    return {
        "peak_overlap_count_per_kb": {
            "available": True,
            "normalization": "peak-overlap count per kb",
        },
        "insulation": {
            "available": insulation_available,
            "normalization": _single_annotation(insulation, "normalization"),
            "balance_status": _single_annotation(
                insulation, "balance_status", "NOT_ASSESSED"
            ),
            "source_resolution_bp": 25_000,
            "reason": None if insulation_available else "No finite insulation values",
        },
        "E1_eigenvector": {
            "available": e1_available,
            "normalization": eigs_status.get("normalization"),
            "balance_status": eigs_status.get("balance_status", "NOT_ASSESSED"),
            "source_resolution_bp": int(eigs_status.get("resolution_bp", 100_000)),
            "reason": (
                None if e1_available
                else eigs_status.get("reason", "No finite E1 values")
            ),
        },
    }


def _bins_for_resolution(clr: cooler.Cooler, chroms: list[str]) -> pd.DataFrame:
    bins = clr.bins()[:]
    bins = bins[bins["chrom"].astype(str).isin(chroms)].reset_index(drop=True)
    bins["bin_idx"] = np.arange(len(bins))
    return bins


def _select_chromosomes(
    chromsizes: dict[str, int],
    drop_chroms: set[str],
    primary_only: bool,
) -> list[str]:
    chroms = [chrom for chrom in chromsizes if chrom not in drop_chroms]
    if primary_only:
        chroms = [chrom for chrom in chroms if PRIMARY_CHROMOSOME.fullmatch(chrom)]
    return chroms


def _peak_overlap_per_bin(bed_path: str | Path, bins: pd.DataFrame) -> np.ndarray:
    """
    Count overlapping MACS3 peak intervals per genomic bin. This is a prototype
    categorical/overlap feature, not continuous ChIP/CUT&Tag signal.
    """
    out = np.zeros(len(bins), dtype=np.float32)
    p = Path(bed_path)
    if not p.exists() or p.stat().st_size == 0:
        return out
    peaks = pd.read_csv(
        p, sep="\t", header=None, usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
    )
    peaks["start"] = pd.to_numeric(peaks["start"], errors="raise")
    peaks["end"] = pd.to_numeric(peaks["end"], errors="raise")
    if (peaks["end"] <= peaks["start"]).any():
        raise ValueError(f"peak BED contains end <= start intervals: {p}")
    for chrom, sub in bins.groupby("chrom", observed=True):
        peaks_c = peaks[peaks["chrom"].astype(str) == str(chrom)]
        if peaks_c.empty:
            continue
        ordered = sub.sort_values("start", kind="stable")
        starts = ordered["start"].to_numpy(dtype=np.int64)
        ends = ordered["end"].to_numpy(dtype=np.int64)
        idx = ordered["bin_idx"].to_numpy(dtype=np.int64)
        if np.any(starts[1:] < ends[:-1]):
            raise ValueError(f"Bins overlap on {chrom}; peak projection requires disjoint bins")

        # Every interval adds one over a contiguous range of disjoint bins. A
        # difference array makes this O(peaks log bins + bins), rather than
        # rescanning every chromosome-wide bin vector once per peak.
        delta = np.zeros(len(ordered) + 1, dtype=np.int64)
        left = np.searchsorted(
            ends, peaks_c["start"].to_numpy(dtype=np.int64), side="right"
        )
        right = np.searchsorted(
            starts, peaks_c["end"].to_numpy(dtype=np.int64), side="left"
        )
        valid = left < right
        np.add.at(delta, left[valid], 1)
        np.add.at(delta, right[valid], -1)
        out[idx] = np.cumsum(delta[:-1], dtype=np.int64).astype(np.float32)
    bin_sizes = (bins["end"].values - bins["start"].values).astype(np.float32)
    return out / np.maximum(bin_sizes / 1000.0, 1.0)


def _interval_feature_per_bin(
    intervals: pd.DataFrame,
    bins: pd.DataFrame,
    value_col: str,
) -> np.ndarray:
    """Project an interval-valued track onto genomic bins.

    Values are averaged by the number of overlapping base pairs. This makes a
    coarse source interval broadcast to *every* finer bin it covers, while a
    coarse destination bin receives the overlap-weighted mean of finer source
    intervals. An exact ``chrom/start`` merge cannot do the former: for
    example, a 100 kb E1 interval would populate only the first of twenty 5 kb
    bins and silently turn the other nineteen into missing values.

    Uncovered bins are returned as NaN so callers can choose their own missing
    value policy.
    """
    required_intervals = {"chrom", "start", "end", value_col}
    required_bins = {"chrom", "start", "end"}
    missing_intervals = required_intervals.difference(intervals.columns)
    missing_bins = required_bins.difference(bins.columns)
    if missing_intervals:
        raise ValueError(f"Interval track is missing required columns: {sorted(missing_intervals)}")
    if missing_bins:
        raise ValueError(f"Bins are missing required columns: {sorted(missing_bins)}")

    out = np.full(len(bins), np.nan, dtype=np.float64)
    if intervals.empty or bins.empty:
        return out.astype(np.float32)

    source = intervals[["chrom", "start", "end", value_col]].copy()
    source["chrom"] = source["chrom"].astype(str)
    for col in ("start", "end", value_col):
        source[col] = pd.to_numeric(source[col], errors="coerce")
    source = source.dropna(subset=["start", "end", value_col])
    source = source[
        np.isfinite(source[["start", "end", value_col]].to_numpy(dtype=float)).all(axis=1)
    ]
    source = source[source["end"] > source["start"]]
    if source.empty:
        return out.astype(np.float32)

    bin_chroms = bins["chrom"].astype(str).to_numpy()
    bin_starts = pd.to_numeric(bins["start"], errors="coerce").to_numpy(dtype=float)
    bin_ends = pd.to_numeric(bins["end"], errors="coerce").to_numpy(dtype=float)
    if (
        not np.isfinite(bin_starts).all()
        or not np.isfinite(bin_ends).all()
        or np.any(bin_ends <= bin_starts)
    ):
        raise ValueError("Bins must have finite coordinates with end > start")

    for chrom in pd.unique(bin_chroms):
        bin_pos = np.flatnonzero(bin_chroms == chrom)
        src = source[source["chrom"] == chrom]
        if not len(bin_pos) or src.empty:
            continue

        bin_order = np.argsort(bin_starts[bin_pos], kind="stable")
        ordered_pos = bin_pos[bin_order]
        bs = bin_starts[ordered_pos]
        be = bin_ends[ordered_pos]

        src = src.sort_values(["start", "end"], kind="stable")
        ss = src["start"].to_numpy(dtype=float)
        se = src["end"].to_numpy(dtype=float)
        sv = src[value_col].to_numpy(dtype=float)
        if np.any(bs[1:] < be[:-1]):
            raise ValueError(f"Bins overlap on {chrom}; interval projection requires disjoint bins")
        if np.any(ss[1:] < se[:-1]):
            raise ValueError(
                f"Source intervals overlap on {chrom}; expected one value per genomic position"
            )

        weighted_sum = np.zeros(len(ordered_pos), dtype=np.float64)
        covered_bp = np.zeros(len(ordered_pos), dtype=np.float64)
        i = j = 0
        while i < len(ordered_pos) and j < len(src):
            overlap = min(be[i], se[j]) - max(bs[i], ss[j])
            if overlap > 0:
                weighted_sum[i] += sv[j] * overlap
                covered_bp[i] += overlap

            # Store both ends before advancing either pointer. A source interval
            # may span many fine bins, and a coarse bin may span many source
            # intervals; equal ends advance both.
            bin_end = be[i]
            source_end = se[j]
            if bin_end <= source_end:
                i += 1
            if source_end <= bin_end:
                j += 1

        has_value = covered_bp > 0
        out[ordered_pos[has_value]] = weighted_sum[has_value] / covered_bp[has_value]

    return out.astype(np.float32)


def _blacklist_mask_for_bins(blacklist_path: str | Path | None, bins: pd.DataFrame) -> np.ndarray:
    """
    Boolean mask (len == n bins) flagging bins that overlap a blacklist interval.

    The blacklist is a required, provenance-tracked input. Missing or malformed
    input therefore fails the export instead of silently claiming that masking was
    applied while returning an all-False mask.
    """
    mask = np.zeros(len(bins), dtype=bool)
    if not blacklist_path or not str(blacklist_path):
        raise ValueError("ORACLE export requires a blacklist path")
    p = Path(blacklist_path)
    if not p.exists() or p.stat().st_size == 0:
        raise FileNotFoundError(f"blacklist is missing or empty: {p}")
    bl = pd.read_csv(
        p, sep="\t", header=None, comment="#", usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
    )
    for col in ("start", "end"):
        bl[col] = pd.to_numeric(bl[col], errors="raise")
    if (bl["end"] <= bl["start"]).any():
        raise ValueError(f"blacklist contains end <= start intervals: {p}")
    for chrom, sub in bins.groupby("chrom", observed=True):
        bl_c = bl[bl["chrom"].astype(str) == str(chrom)]
        if bl_c.empty:
            continue
        starts = sub["start"].values
        ends = sub["end"].values
        idx = sub["bin_idx"].values
        for s, e in zip(bl_c["start"].values, bl_c["end"].values):
            hit = (starts < e) & (ends > s)
            mask[idx[hit]] = True
    return mask


def _as_float(value, default: float) -> float:
    try:
        if pd.isna(value):
            return default
        return float(value)
    except Exception:
        return default


def _loops_to_edges(loops: pd.DataFrame, bins: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """
    Map and consolidate BEDPE loops at one graph resolution.

    Multiple fine-resolution calls can collapse onto the same coarse node pair.
    Exporting them as parallel edges silently weights that pair by call multiplicity
    in generic GNN aggregation. Emit one undirected pair (stored in both directions)
    with max score, min p/q, coarse distance, and the number of collapsed fine calls.
    Min p/q are descriptive strongest-evidence summaries, not re-tested coarse-bin
    significance values.
    """
    if loops.empty:
        return np.zeros((2, 0), dtype=np.int64), np.zeros((0, 5), dtype=np.float32)
    bin_lookup = bins.set_index(["chrom", "start"])["bin_idx"].to_dict()

    grouped: dict[tuple[int, int], list[list[float]]] = {}
    bin_sz = int(bins.iloc[0]["end"] - bins.iloc[0]["start"])
    for _, r in loops.iterrows():
        if str(r["chrom1"]) != str(r["chrom2"]):
            continue
        s1 = int(r["start1"]) // bin_sz * bin_sz
        s2 = int(r["start2"]) // bin_sz * bin_sz
        i = bin_lookup.get((str(r["chrom1"]), s1))
        j = bin_lookup.get((str(r["chrom2"]), s2))
        if i is None or j is None or i == j:
            continue
        key = tuple(sorted((int(i), int(j))))
        grouped.setdefault(key, []).append([
            _as_float(r.get("score", 0.0), 0.0),
            _as_float(r.get("pvalue", 1.0), 1.0),
            _as_float(r.get("fdr", 1.0), 1.0),
            float(abs(s1 - s2)),
        ])

    srcs: list[int] = []
    dsts: list[int] = []
    attrs: list[list[float]] = []
    for (i, j), values in sorted(grouped.items()):
        array = np.asarray(values, dtype=float)
        attr = [
            float(np.max(array[:, 0])),
            float(np.min(array[:, 1])),
            float(np.min(array[:, 2])),
            float(array[0, 3]),
            float(len(values)),
        ]
        srcs.extend([i, j])
        dsts.extend([j, i])
        attrs.extend([attr, attr])
    # A non-empty loop set can contribute no edges when all contacts are trans or
    # collapse into self-edges at a coarse resolution. Preserve a 2-D empty schema.
    if not attrs:
        return np.zeros((2, 0), dtype=np.int64), np.zeros((0, 5), dtype=np.float32)
    return np.array([srcs, dsts], dtype=np.int64), np.asarray(attrs, dtype=np.float32)


def _adjacency_edges(bins: pd.DataFrame) -> np.ndarray:
    """Genomic neighbour edges: bin i ↔ bin i+1 within the same chromosome."""
    idxs: list[tuple[int, int]] = []
    for _, sub in bins.groupby("chrom", observed=True):
        ids = sub["bin_idx"].values
        for a, b in zip(ids[:-1], ids[1:]):
            idxs.append((int(a), int(b)))
            idxs.append((int(b), int(a)))
    if not idxs:
        return np.zeros((2, 0), dtype=np.int64)
    return np.asarray(idxs, dtype=np.int64).T


def _containment_edges(fine_bins: pd.DataFrame, coarse_bins: pd.DataFrame) -> np.ndarray:
    """Map every fine node to exactly one containing coarse node."""
    fine_size = int(fine_bins.iloc[0]["end"] - fine_bins.iloc[0]["start"])
    coarse_size = int(coarse_bins.iloc[0]["end"] - coarse_bins.iloc[0]["start"])
    if coarse_size % fine_size:
        raise ValueError(
            f"coarse resolution {coarse_size} is not divisible by {fine_size}"
        )
    lookup = coarse_bins.set_index(["chrom", "start"])["bin_idx"].to_dict()
    fine_ids = []
    coarse_ids = []
    for row in fine_bins.itertuples(index=False):
        parent_start = int(row.start) // coarse_size * coarse_size
        parent = lookup.get((str(row.chrom), parent_start))
        if parent is None:
            raise ValueError(
                f"no {coarse_size}-bp parent for {row.chrom}:{row.start}-{row.end}"
            )
        fine_ids.append(int(row.bin_idx))
        coarse_ids.append(int(parent))
    return np.asarray([fine_ids, coarse_ids], dtype=np.int64)


def _load_global_tokens(
    microbiome_tsv: str | Path | None,
    sample_id: str,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Load a table-wide, batch-stable optional numeric token schema.

    Columns that contain numeric observations define the key order for every
    sample. Missing cells and absent samples are zero-filled with a false observed
    mask, so graph batching never changes tensor width or loses boundaries.
    """
    if not microbiome_tsv or not str(microbiome_tsv) or not Path(microbiome_tsv).exists():
        return (
            [], np.empty(0, dtype=np.float32), np.empty(0, dtype=bool)
        )
    df = pd.read_csv(
        microbiome_tsv, sep="\t", dtype=str, keep_default_na=False
    )
    if "sample_id" not in df.columns:
        raise ValueError("optional microbiome token table requires a sample_id column")
    if df["sample_id"].duplicated().any():
        raise ValueError("optional microbiome token table has duplicate sample_id values")
    df = df.set_index("sample_id")
    keys: list[str] = []
    numeric: dict[str, pd.Series] = {}
    for column in df.columns:
        source = df[column].astype(str).str.strip()
        values = pd.to_numeric(source.replace("", np.nan), errors="coerce")
        observed_source = source != ""
        if not values.notna().any():
            continue
        if (observed_source & values.isna()).any():
            raise ValueError(
                f"optional token column {column!r} mixes numeric and non-numeric values"
            )
        keys.append(str(column))
        numeric[str(column)] = values

    if not keys:
        return [], np.empty(0, dtype=np.float32), np.empty(0, dtype=bool)
    if sample_id not in df.index:
        return (
            keys, np.zeros(len(keys), dtype=np.float32),
            np.zeros(len(keys), dtype=bool),
        )
    row = np.asarray([numeric[key].loc[sample_id] for key in keys], dtype=float)
    observed = np.isfinite(row)
    values = np.where(observed, row, 0.0).astype(np.float32)
    return keys, values, observed.astype(bool)


def _zero_unobserved_values(x: np.ndarray, observed: np.ndarray) -> np.ndarray:
    """Return a copy with every unavailable feature cell safely zero-filled."""
    if x.shape != observed.shape:
        raise ValueError("feature matrix and observed mask must have identical shape")
    out = np.asarray(x, dtype=np.float32).copy()
    out[~observed.astype(bool)] = 0.0
    return out


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])

    bin_sizes_bp: list[int] = list(snakemake.params.bin_sizes_bp)
    chromsizes = read_chromsizes(snakemake.input.chromsizes)
    drop_chroms = set(snakemake.params.drop_chroms or [])
    primary_only = bool(snakemake.params.primary_chromosomes_only)
    chroms = _select_chromosomes(chromsizes, drop_chroms, primary_only)
    if not chroms:
        raise ValueError("ORACLE chromosome policy selected no chromosomes")

    loops_annot_path = Path(snakemake.input.loops_annot)
    loops = load_loops_bedpe(loops_annot_path)
    peaks = Path(snakemake.input.peaks)

    insul = pd.read_csv(snakemake.input.insul, sep="\t")
    eigs = pd.read_csv(snakemake.input.eigs, sep="\t")
    eigs_status = json.loads(Path(snakemake.input.eigs_status).read_text())
    balance = json.loads(Path(snakemake.input.balance).read_text())
    ix_col = select_insulation_column(insul)
    feature_availability = node_feature_availability(
        insul, ix_col, eigs, eigs_status
    )
    balance_sha256 = _sha256(snakemake.input.balance)

    sample_id = snakemake.wildcards.sample
    assembly = str(snakemake.params.assembly)
    mark = str(snakemake.params.mark)
    cell_type = str(snakemake.params.cell_type)
    blacklist_path = str(snakemake.input.blacklist)
    microbiome_inputs = list(snakemake.input.microbiome)
    microbiome_tsv = microbiome_inputs[0] if microbiome_inputs else None
    micro_keys, micro_values, micro_observed = _load_global_tokens(
        microbiome_tsv, sample_id
    )
    n_blacklisted_bins: dict[str, int] = {}

    Path(snakemake.output.h5).parent.mkdir(parents=True, exist_ok=True)
    graphs: dict[str, dict] = {}
    h5 = h5py.File(snakemake.output.h5, "w")
    h5.attrs["schema"] = "oracle-hichip-cos-v2"
    h5.attrs["sample_id"] = sample_id
    h5.attrs["assembly"] = assembly
    h5.attrs["mark"] = mark
    h5.attrs["cell_type"] = cell_type
    h5.attrs["coordinate_system"] = "0-based half-open"
    h5.attrs["node_feature_channels"] = json.dumps(NODE_FEATURE_CHANNELS)
    h5.attrs["node_feature_availability"] = json.dumps(feature_availability)
    h5.attrs["balance_qc_status"] = str(balance.get("status", "NOT_ASSESSED"))
    h5.attrs["balance_qc_sha256"] = balance_sha256
    h5.attrs["loop_edge_attr_channels"] = json.dumps(
        [
            "max_loop_score", "min_loop_pvalue", "min_loop_fdr",
            "genomic_distance_bp", "fine_loop_count",
        ]
    )

    for bp in bin_sizes_bp:
        clr = cooler.Cooler(f"{snakemake.input.mcool}::resolutions/{bp}")
        bins = _bins_for_resolution(clr, chroms)

        peak_overlap_count = _peak_overlap_per_bin(peaks, bins)

        # Zero the peak-overlap feature for bins overlapping ENCODE blacklist
        # regions — these produce mapping artefacts and spurious peak overlaps.
        bl_mask = _blacklist_mask_for_bins(blacklist_path, bins)
        if bl_mask.any():
            peak_overlap_count[bl_mask] = 0.0
        n_blacklisted_bins[f"res_{bp}"] = int(bl_mask.sum())

        insulation_raw = _interval_feature_per_bin(insul, bins, ix_col)
        e1_raw = _interval_feature_per_bin(eigs, bins, "E1")
        insulation = np.nan_to_num(insulation_raw, nan=0.0).astype(np.float32)
        e1 = np.nan_to_num(e1_raw, nan=0.0).astype(np.float32)

        # Prototype node features only; additional modalities require a documented,
        # coordinate-compatible integration step.
        x = np.stack([peak_overlap_count, insulation, e1], axis=1).astype(np.float32)
        x_observed_mask = np.stack(
            [
                ~bl_mask,
                np.isfinite(insulation_raw) & ~bl_mask,
                np.isfinite(e1_raw)
                & bool(feature_availability["E1_eigenvector"]["available"])
                & ~bl_mask,
            ],
            axis=1,
        ).astype(bool)
        # A numeric zero plus an explicit false mask is safer than retaining a
        # plausible-looking value in a blacklisted or unavailable feature cell.
        x = _zero_unobserved_values(x, x_observed_mask)

        edge_loops, edge_attr_loops = _loops_to_edges(loops, bins)
        edge_adj = _adjacency_edges(bins)

        grp = h5.create_group(f"res_{bp}")
        grp.create_dataset("bin_chrom", data=bins["chrom"].astype(str).values.astype("S"))
        grp.create_dataset("bin_start", data=bins["start"].values.astype(np.int64))
        grp.create_dataset("bin_end", data=bins["end"].values.astype(np.int64))
        grp.create_dataset("x", data=x, compression="gzip", compression_opts=4)
        grp.create_dataset(
            "x_observed_mask", data=x_observed_mask,
            compression="gzip", compression_opts=4,
        )
        grp.create_dataset("blacklist_mask", data=bl_mask.astype(bool))
        grp.create_dataset("loop_edge_index", data=edge_loops, compression="gzip", compression_opts=4)
        grp.create_dataset("loop_edge_attr", data=edge_attr_loops, compression="gzip", compression_opts=4)
        grp.create_dataset("adjacency_edge_index", data=edge_adj, compression="gzip", compression_opts=4)
        grp.attrs["n_nodes"] = x.shape[0]
        grp.attrs["n_loop_edges"] = edge_loops.shape[1]
        grp.attrs["n_adjacency_edges"] = edge_adj.shape[1]
        grp.attrs["bin_size_bp"] = bp

        chrom_names = list(pd.unique(bins["chrom"].astype(str)))
        chrom_to_id = {chrom: i for i, chrom in enumerate(chrom_names)}
        graphs[f"res_{bp}"] = {
            "x": x,
            "x_observed_mask": x_observed_mask,
            "blacklist_mask": bl_mask.astype(bool),
            "loop_edge_index": edge_loops,
            "loop_edge_attr": edge_attr_loops,
            "adjacency_edge_index": edge_adj,
            "bins": bins,
            "chrom_id": bins["chrom"].astype(str).map(chrom_to_id).to_numpy(dtype=np.int16),
            "chrom_names": chrom_names,
            "bin_start": bins["start"].to_numpy(dtype=np.int64),
            "bin_end": bins["end"].to_numpy(dtype=np.int64),
        }

    hierarchy_edges: dict[tuple[int, int], np.ndarray] = {}
    hierarchy = h5.create_group("hierarchy")
    ordered_resolutions = sorted(bin_sizes_bp)
    for fine_bp, coarse_bp in zip(ordered_resolutions[:-1], ordered_resolutions[1:]):
        edges = _containment_edges(
            graphs[f"res_{fine_bp}"]["bins"],
            graphs[f"res_{coarse_bp}"]["bins"],
        )
        hierarchy_edges[(fine_bp, coarse_bp)] = edges
        relation = hierarchy.create_group(f"res_{fine_bp}_to_res_{coarse_bp}")
        relation.create_dataset(
            "contained_by_edge_index", data=edges,
            compression="gzip", compression_opts=4,
        )
        relation.create_dataset(
            "contains_edge_index", data=edges[[1, 0], :],
            compression="gzip", compression_opts=4,
        )

    h5.attrs["microbiome_keys"] = json.dumps(micro_keys)
    h5.create_dataset("microbiome_values", data=micro_values.reshape(1, -1))
    h5.create_dataset(
        "microbiome_observed_mask", data=micro_observed.reshape(1, -1)
    )
    h5.close()

    try:
        import torch
        from torch_geometric.data import HeteroData
    except Exception as exc:
        raise RuntimeError(
            "PyTorch Geometric is required to create the ORACLE .pt output. "
            "Install torch and torch_geometric, or change the Snakemake outputs "
            "before running an HDF5-only export."
        ) from exc

    data = HeteroData()
    data.schema = "oracle-hichip-cos-v2"
    data.coordinate_system = "0-based half-open"
    data.node_feature_channels = NODE_FEATURE_CHANNELS
    data.node_feature_availability = feature_availability
    data.loop_edge_attr_channels = [
        "max_loop_score", "min_loop_pvalue", "min_loop_fdr",
        "genomic_distance_bp", "fine_loop_count",
    ]
    data.balance_qc_status = str(balance.get("status", "NOT_ASSESSED"))
    data.balance_qc_sha256 = balance_sha256
    data["sample"].id = sample_id
    data["sample"].assembly = assembly
    data["sample"].mark = mark
    data["sample"].cell_type = cell_type
    data["sample"].num_nodes = 1
    data["sample"].microbiome = torch.from_numpy(
        micro_values.reshape(1, -1)
    ).float()
    data["sample"].microbiome_observed_mask = torch.from_numpy(
        micro_observed.reshape(1, -1)
    ).bool()

    for key, g in graphs.items():
        node_type = f"bin_{key}"
        data[node_type].x = torch.from_numpy(g["x"])
        data[node_type].x_observed_mask = torch.from_numpy(
            g["x_observed_mask"]
        ).bool()
        data[node_type].blacklist_mask = torch.from_numpy(
            g["blacklist_mask"]
        ).bool()
        data[node_type].chrom_id = torch.from_numpy(g["chrom_id"]).long()
        data[node_type].chrom_names = g["chrom_names"]
        data[node_type].bin_start = torch.from_numpy(g["bin_start"]).long()
        data[node_type].bin_end = torch.from_numpy(g["bin_end"]).long()
        loop_relation = (node_type, "loop", node_type)
        data[loop_relation].edge_index = torch.from_numpy(
            g["loop_edge_index"]
        ).long()
        data[loop_relation].edge_attr = torch.from_numpy(g["loop_edge_attr"])
        adjacency_relation = (node_type, "adjacent", node_type)
        data[adjacency_relation].edge_index = torch.from_numpy(
            g["adjacency_edge_index"]
        ).long()

    for (fine_bp, coarse_bp), edges in hierarchy_edges.items():
        fine_type = f"bin_res_{fine_bp}"
        coarse_type = f"bin_res_{coarse_bp}"
        data[(fine_type, "contained_by", coarse_type)].edge_index = (
            torch.from_numpy(edges).long()
        )
        data[(coarse_type, "contains", fine_type)].edge_index = (
            torch.from_numpy(edges[[1, 0], :]).long()
        )

    torch.save(data, snakemake.output.pt)

    input_hashes = {
        "mcool":        _sha256(snakemake.input.mcool),
        "loops_annot":  _sha256(snakemake.input.loops_annot),
        "peaks":        _sha256(snakemake.input.peaks),
        "insulation":   _sha256(snakemake.input.insul),
        "eigenvectors": _sha256(snakemake.input.eigs),
        "eigenvector_status": _sha256(snakemake.input.eigs_status),
        "balance_qc":    balance_sha256,
        "loop_qc":      _sha256(snakemake.input.loop_qc),
        "chromsizes":    _sha256(snakemake.input.chromsizes),
        "blacklist":     _sha256(snakemake.input.blacklist),
    }
    for index, path in enumerate(snakemake.input.shared_code):
        input_hashes[f"shared_code_{index}_{Path(path).name}"] = _sha256(path)
    if microbiome_tsv:
        input_hashes["microbiome_tokens"] = _sha256(microbiome_tsv)

    blacklist_applied = True
    manifest = {
        "schema": "oracle-hichip-cos-v2",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sample_id": sample_id,
        "assembly": assembly,
        "mark": mark,
        "cell_type": cell_type,
        "coordinate_system": "0-based half-open",
        "node_order": "resolution-local cooler bin order after configured chromosome filtering",
        "chromosome_policy": (
            "primary chr1-22/X only" if primary_only
            else "all chromsizes entries except configured exclusions"
        ),
        "included_chromosomes": chroms,
        "resolutions_bp": bin_sizes_bp,
        "node_feature_channels": NODE_FEATURE_CHANNELS,
        "node_feature_availability": feature_availability,
        "unavailable_feature_value_policy": (
            "Numeric channel order is fixed; unavailable source values are zero-filled. "
            "Consumers must apply each node store's x_observed_mask and blacklist_mask "
            "and consult node_feature_availability before use."
        ),
        "balance_qc": {
            "path": str(snakemake.input.balance),
            "sha256": balance_sha256,
            "status": balance.get("status", "NOT_ASSESSED"),
            "n_configured": balance.get("n_configured"),
            "n_converged": balance.get("n_converged"),
            "n_nonconverged": balance.get("n_nonconverged"),
            "n_missing": balance.get("n_missing"),
            "nonconverged_resolutions_bp": balance.get(
                "nonconverged_resolutions_bp", []
            ),
            "missing_resolutions_bp": balance.get("missing_resolutions_bp", []),
        },
        "relations": {
            "within_resolution_edge_types": [
                [f"bin_res_{bp}", relation, f"bin_res_{bp}"]
                for bp in sorted(bin_sizes_bp)
                for relation in ("loop", "adjacent")
            ],
            "hierarchy_edge_types": [
                edge_type
                for fine_bp, coarse_bp in zip(
                    sorted(bin_sizes_bp)[:-1], sorted(bin_sizes_bp)[1:]
                )
                for edge_type in (
                    [f"bin_res_{fine_bp}", "contained_by", f"bin_res_{coarse_bp}"],
                    [f"bin_res_{coarse_bp}", "contains", f"bin_res_{fine_bp}"],
                )
            ],
        },
        "loop_edge_attr_channels": [
            "max_loop_score", "min_loop_pvalue", "min_loop_fdr",
            "genomic_distance_bp", "fine_loop_count",
        ],
        "coarse_loop_aggregation": (
            "One edge per unordered coarse node pair, stored in both directions. "
            "Score=max, p/q=min, and fine_loop_count records collapsed calls; "
            "coarse p/q are descriptive summaries, not re-tested significance."
        ),
        "hierarchy": (
            "Every fine node has one contained_by edge to the next configured "
            "coarser resolution and a reverse contains edge."
        ),
        "blacklist": str(blacklist_path),
        "blacklist_applied": blacklist_applied,
        "blacklist_sha256": _sha256(str(blacklist_path)),
        "n_blacklisted_bins_per_resolution": n_blacklisted_bins,
        "microbiome_keys": micro_keys,
        "microbiome_observed_mask": micro_observed.tolist(),
        "scope_note": (
            "Prototype HiChIP COS export: peak-overlap count is not continuous per-mark "
            "ChIP/CUT&Tag signal. Add modality channels only through a documented, "
            "coordinate-compatible integration step."
        ),
        "input_sha256": input_hashes,
        "outputs": {
            "pt": {"path": str(snakemake.output.pt), "sha256": _sha256(snakemake.output.pt)},
            "h5": {"path": str(snakemake.output.h5), "sha256": _sha256(snakemake.output.h5)},
        },
    }
    Path(snakemake.output.manifest).write_text(json.dumps(manifest, indent=2))


# Guarded so the module can be imported by the tests. Snakemake injects
# `snakemake` into the script's globals before executing it.
if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
