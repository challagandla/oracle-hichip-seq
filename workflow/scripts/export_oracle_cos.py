"""
Export a HiChIP sample to the canonical ORACLE Chromatin Operating System
(COS) representation expected by downstream foundation-model code.

Current scope: HiChIP structural graph + peak-overlap prototype node features.
True continuous per-mark signal tracks (bigWig/RPKM), ATAC and methylation
channels are intentionally left for sister modality pipelines and later merge.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_loops_bedpe, read_chromsizes, setup_logging  # noqa: E402


def _bins_for_resolution(clr: cooler.Cooler, chroms: list[str]) -> pd.DataFrame:
    bins = clr.bins()[:]
    bins = bins[bins["chrom"].astype(str).isin(chroms)].reset_index(drop=True)
    bins["bin_idx"] = np.arange(len(bins))
    return bins


def _peak_overlap_per_bin(bed_path: str | Path, bins: pd.DataFrame) -> np.ndarray:
    """
    Count overlapping MACS2 peak intervals per genomic bin. This is a prototype
    categorical/overlap feature, not continuous ChIP/CUT&Tag signal.
    """
    out = np.zeros(len(bins), dtype=np.float32)
    p = Path(bed_path)
    if not p.exists() or p.stat().st_size == 0:
        return out
    peaks = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"])
    for chrom, sub in bins.groupby("chrom"):
        peaks_c = peaks[peaks["chrom"].astype(str) == str(chrom)]
        if peaks_c.empty:
            continue
        starts = sub["start"].values
        ends = sub["end"].values
        idx = sub["bin_idx"].values
        for s, e in zip(peaks_c["start"].values, peaks_c["end"].values):
            mask = (starts < e) & (ends > s)
            out[idx[mask]] += 1.0
    bin_sizes = (bins["end"].values - bins["start"].values).astype(np.float32)
    return out / np.maximum(bin_sizes / 1000.0, 1.0)


def _as_float(value, default: float) -> float:
    try:
        if pd.isna(value):
            return default
        return float(value)
    except Exception:
        return default


def _loops_to_edges(loops: pd.DataFrame, bins: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """
    Map BEDPE loops to (edge_index, edge_attr) by bin coordinates. Edges are
    intra-chromosomal only; interchromosomal contacts are excluded from this
    first COS graph representation.
    """
    if loops.empty:
        return np.zeros((2, 0), dtype=np.int64), np.zeros((0, 3), dtype=np.float32)
    bin_lookup = bins.set_index(["chrom", "start"])["bin_idx"].to_dict()

    srcs: list[int] = []
    dsts: list[int] = []
    attrs: list[list[float]] = []
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
        attr = [
            _as_float(r.get("score", 0.0), 0.0),
            _as_float(r.get("fdr", 1.0), 1.0),
            float(abs(s1 - s2)),
        ]
        srcs.extend([i, j])
        dsts.extend([j, i])
        attrs.extend([attr, attr])
    return np.array([srcs, dsts], dtype=np.int64), np.asarray(attrs, dtype=np.float32)


def _adjacency_edges(bins: pd.DataFrame) -> np.ndarray:
    """Genomic neighbour edges: bin i ↔ bin i+1 within the same chromosome."""
    idxs: list[tuple[int, int]] = []
    for _, sub in bins.groupby("chrom"):
        ids = sub["bin_idx"].values
        for a, b in zip(ids[:-1], ids[1:]):
            idxs.append((int(a), int(b)))
            idxs.append((int(b), int(a)))
    if not idxs:
        return np.zeros((2, 0), dtype=np.int64)
    return np.asarray(idxs, dtype=np.int64).T


def _load_global_tokens(microbiome_tsv: str | Path | None, sample_id: str) -> tuple[list[str], np.ndarray]:
    """Optional numeric microbiome/metabolite tokens. Non-numeric columns are skipped."""
    if not microbiome_tsv or not str(microbiome_tsv) or not Path(microbiome_tsv).exists():
        return [], np.empty(0, dtype=np.float32)
    df = pd.read_csv(microbiome_tsv, sep="\t").set_index("sample_id")
    if sample_id not in df.index:
        return [], np.empty(0, dtype=np.float32)
    row = pd.to_numeric(df.loc[sample_id], errors="coerce").dropna()
    return row.index.astype(str).tolist(), row.to_numpy(dtype=np.float32)


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])

    bin_sizes_bp: list[int] = list(snakemake.params.bin_sizes_bp)
    chromsizes = read_chromsizes(snakemake.params.chromsizes)
    drop_chroms = set(snakemake.params.drop_chroms or [])
    chroms = [c for c in chromsizes if c not in drop_chroms]

    loops_annot_path = Path(snakemake.input.loops_annot)
    loops = load_loops_bedpe(loops_annot_path)
    peaks = Path(snakemake.input.peaks)

    insul = pd.read_csv(snakemake.input.insul, sep="\t")
    eigs = pd.read_csv(snakemake.input.eigs, sep="\t")

    sample_id = snakemake.wildcards.sample
    micro_keys, micro_values = _load_global_tokens(snakemake.params.microbiome_tsv, sample_id)

    Path(snakemake.output.h5).parent.mkdir(parents=True, exist_ok=True)
    graphs: dict[str, dict] = {}
    h5 = h5py.File(snakemake.output.h5, "w")
    h5.attrs["sample_id"] = sample_id

    for bp in bin_sizes_bp:
        clr = cooler.Cooler(f"{snakemake.input.mcool}::resolutions/{bp}")
        bins = _bins_for_resolution(clr, chroms)

        peak_overlap_count = _peak_overlap_per_bin(peaks, bins)

        ix = insul.copy()
        ix["start"] = (ix["start"] // bp) * bp
        ix_col = "log2_insulation_score" if "log2_insulation_score" in ix.columns else ix.columns[-1]
        ix = ix.groupby(["chrom", "start"], as_index=False)[ix_col].mean()
        eigs_b = eigs.copy()
        eigs_b["start"] = (eigs_b["start"] // bp) * bp
        eigs_b = eigs_b.groupby(["chrom", "start"], as_index=False)["E1"].mean()

        bins = bins.merge(ix, on=["chrom", "start"], how="left").merge(eigs_b, on=["chrom", "start"], how="left")
        insulation = bins[ix_col].fillna(0.0).to_numpy(dtype=np.float32)
        e1 = bins["E1"].fillna(0.0).to_numpy(dtype=np.float32)

        # Prototype node features only. Sister pipelines should add true signal channels.
        x = np.stack([peak_overlap_count, insulation, e1], axis=1).astype(np.float32)

        edge_loops, edge_attr_loops = _loops_to_edges(loops, bins)
        edge_adj = _adjacency_edges(bins)
        edge_index = np.concatenate([edge_loops, edge_adj], axis=1) if edge_loops.size else edge_adj
        edge_attr = np.concatenate([
            edge_attr_loops,
            np.zeros((edge_adj.shape[1], edge_attr_loops.shape[1] if edge_attr_loops.size else 3), dtype=np.float32),
        ], axis=0)
        edge_type = np.concatenate([
            np.ones(edge_loops.shape[1], dtype=np.int8),
            np.zeros(edge_adj.shape[1], dtype=np.int8),
        ])

        grp = h5.create_group(f"res_{bp}")
        grp.create_dataset("bin_chrom", data=bins["chrom"].astype(str).values.astype("S"))
        grp.create_dataset("bin_start", data=bins["start"].values.astype(np.int64))
        grp.create_dataset("bin_end", data=bins["end"].values.astype(np.int64))
        grp.create_dataset("x", data=x, compression="gzip", compression_opts=4)
        grp.create_dataset("edge_index", data=edge_index, compression="gzip", compression_opts=4)
        grp.create_dataset("edge_attr", data=edge_attr, compression="gzip", compression_opts=4)
        grp.create_dataset("edge_type", data=edge_type)
        grp.attrs["n_nodes"] = x.shape[0]
        grp.attrs["n_edges"] = edge_index.shape[1]
        grp.attrs["bin_size_bp"] = bp

        graphs[f"res_{bp}"] = {"x": x, "edge_index": edge_index, "edge_attr": edge_attr, "edge_type": edge_type}

    h5.attrs["microbiome_keys"] = json.dumps(micro_keys)
    if len(micro_values):
        h5.create_dataset("microbiome_values", data=micro_values)
    h5.close()

    try:
        import torch
        from torch_geometric.data import HeteroData
    except Exception as exc:
        Path(snakemake.output.pt).write_bytes(b"")
        Path(snakemake.output.manifest).write_text(json.dumps({"sample_id": sample_id, "torch_geometric": False, "error": str(exc)}, indent=2))
        return

    data = HeteroData()
    data["sample"].id = sample_id
    data["sample"].microbiome = torch.tensor(micro_values, dtype=torch.float32) if len(micro_values) else torch.empty(0)

    for key, g in graphs.items():
        node_type = f"bin_{key}"
        data[node_type].x = torch.from_numpy(g["x"])
        edge_type = (node_type, "contact", node_type)
        data[edge_type].edge_index = torch.from_numpy(g["edge_index"]).long()
        data[edge_type].edge_attr = torch.from_numpy(g["edge_attr"])
        data[edge_type].edge_kind = torch.from_numpy(g["edge_type"])

    torch.save(data, snakemake.output.pt)

    def _sha256(path: str) -> str:
        """Return hex SHA-256 of a file, or 'MISSING' if it does not exist."""
        p = Path(path)
        if not p.exists():
            return "MISSING"
        h = hashlib.sha256()
        with open(p, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()

    input_hashes = {
        "mcool":        _sha256(snakemake.input.mcool),
        "loops_annot":  _sha256(snakemake.input.loops_annot),
        "peaks":        _sha256(snakemake.input.peaks),
        "insulation":   _sha256(snakemake.input.insul),
        "eigenvectors": _sha256(snakemake.input.eigs),
        "loop_qc":      _sha256(snakemake.input.loop_qc),
    }

    manifest = {
        "sample_id": sample_id,
        "resolutions_bp": bin_sizes_bp,
        "node_feature_channels": ["peak_overlap_count_per_kb", "insulation", "E1_eigenvector"],
        "edge_kinds": {"0": "adjacency", "1": "loop"},
        "edge_attr_channels": ["loop_score", "loop_fdr", "genomic_distance_bp"],
        "microbiome_keys": micro_keys,
        "scope_note": (
            "Prototype HiChIP COS export: peak-overlap count is not continuous per-mark "
            "ChIP/CUT&Tag signal. Additional modality channels should be merged by sister pipelines."
        ),
        "input_sha256": input_hashes,
        "outputs": {"pt": str(snakemake.output.pt), "h5": str(snakemake.output.h5)},
    }
    Path(snakemake.output.manifest).write_text(json.dumps(manifest, indent=2))


main(snakemake)  # type: ignore[name-defined]  # noqa: F821
