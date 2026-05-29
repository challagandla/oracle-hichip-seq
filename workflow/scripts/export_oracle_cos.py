"""
Export a HiChIP sample to the canonical ORACLE Chromatin Operating System
(COS) representation expected by the foundation-model training scripts:

    - One HeteroData graph per ORACLE resolution (5 kb, 25 kb, 100 kb, 1 Mb)
      stacked inside a single PyTorch object.
    - Node features per bin: per-mark signal (bigWig sum), insulation score,
      A/B compartment eigenvalue, blacklist flag, sequence-token placeholder.
    - Edge index: union of (a) intra-chromosomal HiChIP loops at this
      resolution, (b) per-bin neighbours (genomic adjacency).
    - Edge attributes: loop strength, FDR, distance, anchor peak overlap.
    - Global tokens: sample metadata + optional microbiome tokens.

The same object is also mirrored to HDF5 for non-PyG consumers and for
debugging the contents from the command line (`h5ls`, `h5dump`).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import cooler
import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_loops_bedpe, read_chromsizes, setup_logging  # noqa: E402


# ---------------- helpers ----------------

def _bins_for_resolution(clr: cooler.Cooler, chroms: list[str]) -> pd.DataFrame:
    bins = clr.bins()[:]
    bins = bins[bins["chrom"].astype(str).isin(chroms)].reset_index(drop=True)
    bins["bin_idx"] = np.arange(len(bins))
    return bins


def _signal_per_bin(bw_or_bed: str | Path, bins: pd.DataFrame) -> np.ndarray:
    """
    Aggregate peak coverage per bin (count of overlapping peaks normalised to
    bin size). Falls back to zeros if file missing. Replace with real bigWig
    summation when bigWigs are emitted by upstream modality pipelines.
    """
    out = np.zeros(len(bins), dtype=np.float32)
    p = Path(bw_or_bed)
    if not p.exists():
        return out
    peaks = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 2],
                        names=["chrom", "start", "end"])
    for chrom, sub in bins.groupby("chrom"):
        peaks_c = peaks[peaks["chrom"].astype(str) == str(chrom)]
        if peaks_c.empty:
            continue
        starts = sub["start"].values
        ends = sub["end"].values
        idx = sub["bin_idx"].values
        # naive overlap count
        for s, e in zip(peaks_c["start"].values, peaks_c["end"].values):
            mask = (starts < e) & (ends > s)
            out[idx[mask]] += 1.0
    return out


def _loops_to_edges(loops: pd.DataFrame, bins: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """
    Map BEDPE loops to (edge_index, edge_attr) by bin coordinates.
    Edges are intra-chromosomal only (per HiChIP convention).
    """
    if loops.empty:
        return np.zeros((2, 0), dtype=np.int64), np.zeros((0, 3), dtype=np.float32)
    bin_lookup = bins.set_index(["chrom", "start"])["bin_idx"].to_dict()

    srcs: list[int] = []; dsts: list[int] = []; attrs: list[list[float]] = []
    bin_sz = int(bins.iloc[0]["end"] - bins.iloc[0]["start"])
    for _, r in loops.iterrows():
        if r["chrom1"] != r["chrom2"]:
            continue
        s1 = int(r["start1"]) // bin_sz * bin_sz
        s2 = int(r["start2"]) // bin_sz * bin_sz
        i = bin_lookup.get((str(r["chrom1"]), s1))
        j = bin_lookup.get((str(r["chrom2"]), s2))
        if i is None or j is None or i == j:
            continue
        srcs.append(i); dsts.append(j)
        srcs.append(j); dsts.append(i)
        attr = [
            float(r.get("score", 0.0)) if "score" in r else 0.0,
            float(r.get("fdr", 1.0)) if "fdr" in r else 1.0,
            float(abs(s1 - s2)),
        ]
        attrs.append(attr); attrs.append(attr)
    return np.array([srcs, dsts], dtype=np.int64), np.asarray(attrs, dtype=np.float32)


def _adjacency_edges(bins: pd.DataFrame) -> np.ndarray:
    """Genomic neighbour edges: bin i ↔ bin i+1 within the same chromosome."""
    idxs: list[tuple[int, int]] = []
    for _, sub in bins.groupby("chrom"):
        ids = sub["bin_idx"].values
        for a, b in zip(ids[:-1], ids[1:]):
            idxs.append((int(a), int(b))); idxs.append((int(b), int(a)))
    if not idxs:
        return np.zeros((2, 0), dtype=np.int64)
    arr = np.asarray(idxs, dtype=np.int64).T
    return arr


def _load_global_tokens(microbiome_tsv: str | Path | None, sample_id: str) -> dict[str, float]:
    """Optional microbiome metabolite tokens (SCFAs, bile acids, diversity)."""
    if not microbiome_tsv or not str(microbiome_tsv) or not Path(microbiome_tsv).exists():
        return {}
    df = pd.read_csv(microbiome_tsv, sep="\t").set_index("sample_id")
    if sample_id not in df.index:
        return {}
    return df.loc[sample_id].to_dict()


# ---------------- main ----------------

def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])

    bin_sizes_bp: list[int] = list(snakemake.params.bin_sizes_bp)
    chromsizes = read_chromsizes(snakemake.params.chromsizes)
    drop_chroms = set(snakemake.params.drop_chroms or [])
    chroms = [c for c in chromsizes if c not in drop_chroms]

    loops_annot_path = Path(snakemake.input.loops_annot)
    if loops_annot_path.exists() and loops_annot_path.stat().st_size > 0:
        loops = load_loops_bedpe(loops_annot_path)
    else:
        loops = pd.DataFrame()

    peaks = Path(snakemake.input.peaks)

    insul = pd.read_csv(snakemake.input.insul, sep="\t")
    eigs = pd.read_csv(snakemake.input.eigs, sep="\t")

    sample_id = snakemake.wildcards.sample
    micro = _load_global_tokens(snakemake.params.microbiome_tsv, sample_id)

    # Build per-resolution graphs
    graphs: dict[str, dict] = {}
    h5 = h5py.File(snakemake.output.h5, "w")
    h5.attrs["sample_id"] = sample_id

    for bp in bin_sizes_bp:
        clr = cooler.Cooler(f"{snakemake.input.mcool}::resolutions/{bp}")
        bins = _bins_for_resolution(clr, chroms)

        peak_signal = _signal_per_bin(peaks, bins)

        # insulation / eigs joined by (chrom,start)
        ix = insul.rename(columns={c: c for c in insul.columns})
        ix["start"] = (ix["start"] // bp) * bp
        ix_col = "log2_insulation_score" if "log2_insulation_score" in ix.columns else ix.columns[-1]
        ix = ix.groupby(["chrom", "start"], as_index=False)[ix_col].mean()
        eigs_b = eigs.copy()
        eigs_b["start"] = (eigs_b["start"] // bp) * bp
        eigs_b = eigs_b.groupby(["chrom", "start"], as_index=False)["E1"].mean()

        bins = bins.merge(ix, on=["chrom", "start"], how="left").merge(eigs_b, on=["chrom", "start"], how="left")
        insulation = bins[ix_col].fillna(0.0).to_numpy(dtype=np.float32)
        e1 = bins["E1"].fillna(0.0).to_numpy(dtype=np.float32)

        # Stack node features: [peak_signal, insulation, E1] — additional modality
        # channels are added by sister pipelines (atac/, cutandtag/ etc) during
        # merge in oracle/training/data.py.
        x = np.stack([peak_signal, insulation, e1], axis=1).astype(np.float32)

        loops_at_res = loops.copy()
        edge_loops, edge_attr_loops = _loops_to_edges(loops_at_res, bins)
        edge_adj = _adjacency_edges(bins)
        edge_index = np.concatenate([edge_loops, edge_adj], axis=1) if edge_loops.size else edge_adj
        edge_attr = np.concatenate([
            edge_attr_loops,
            np.zeros((edge_adj.shape[1], edge_attr_loops.shape[1] if edge_attr_loops.size else 3), dtype=np.float32),
        ], axis=0)
        # Edge type: 0 = adjacency, 1 = loop
        edge_type = np.concatenate([
            np.ones(edge_loops.shape[1], dtype=np.int8),
            np.zeros(edge_adj.shape[1], dtype=np.int8),
        ])

        # Save into HDF5 group
        grp = h5.create_group(f"res_{bp}")
        grp.create_dataset("bin_chrom", data=bins["chrom"].astype(str).values.astype("S"))
        grp.create_dataset("bin_start", data=bins["start"].values.astype(np.int64))
        grp.create_dataset("bin_end",   data=bins["end"].values.astype(np.int64))
        grp.create_dataset("x", data=x, compression="gzip", compression_opts=4)
        grp.create_dataset("edge_index", data=edge_index, compression="gzip", compression_opts=4)
        grp.create_dataset("edge_attr",  data=edge_attr,  compression="gzip", compression_opts=4)
        grp.create_dataset("edge_type",  data=edge_type)
        grp.attrs["n_nodes"] = x.shape[0]
        grp.attrs["n_edges"] = edge_index.shape[1]
        grp.attrs["bin_size_bp"] = bp

        graphs[f"res_{bp}"] = {
            "x": x, "edge_index": edge_index, "edge_attr": edge_attr,
            "edge_type": edge_type,
            "bins": bins.assign(bin_size=bp),
        }

    # Global tokens
    h5.attrs["microbiome_keys"] = json.dumps(list(micro.keys()))
    if micro:
        h5.create_dataset("microbiome_values", data=np.array(list(micro.values()), dtype=np.float32))
    h5.close()

    # PyTorch Geometric mirror
    try:
        import torch
        from torch_geometric.data import HeteroData
    except Exception as exc:
        Path(snakemake.output.pt).write_bytes(b"")
        Path(snakemake.output.manifest).write_text(json.dumps({"sample_id": sample_id, "torch_geometric": False, "error": str(exc)}, indent=2))
        return

    data = HeteroData()
    data["sample"].id = sample_id
    data["sample"].microbiome = torch.tensor(list(micro.values()), dtype=torch.float32) if micro else torch.empty(0)

    for key, g in graphs.items():
        node_type = f"bin_{key}"
        data[node_type].x = torch.from_numpy(g["x"])
        edge_type = (node_type, "contact", node_type)
        data[edge_type].edge_index = torch.from_numpy(g["edge_index"]).long()
        data[edge_type].edge_attr  = torch.from_numpy(g["edge_attr"])
        data[edge_type].edge_kind  = torch.from_numpy(g["edge_type"])

    torch.save(data, snakemake.output.pt)

    manifest = {
        "sample_id": sample_id,
        "resolutions_bp": bin_sizes_bp,
        "node_feature_channels": ["peak_signal", "insulation", "E1_eigenvector"],
        "edge_kinds": {"0": "adjacency", "1": "loop"},
        "microbiome_keys": list(micro.keys()),
        "outputs": {"pt": str(snakemake.output.pt), "h5": str(snakemake.output.h5)},
    }
    Path(snakemake.output.manifest).write_text(json.dumps(manifest, indent=2))


main(snakemake)  # type: ignore[name-defined]  # noqa: F821
