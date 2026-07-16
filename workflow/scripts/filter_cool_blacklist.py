"""Remove every contact touching a blacklist-overlapping bin from a cooler."""
import hashlib
import os
import sys
import uuid
from pathlib import Path

import cooler
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from utils import open_text_auto, setup_logging, write_json  # noqa: E402


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_blacklist(path: str | Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    intervals: dict[str, list[tuple[int, int]]] = {}
    with open_text_auto(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end = line.rstrip().split("\t")[:3]
            intervals.setdefault(chrom, []).append((int(start), int(end)))
    index = {}
    for chrom, values in intervals.items():
        ordered = np.asarray(sorted(values), dtype=np.int64)
        index[chrom] = (ordered[:, 0], np.maximum.accumulate(ordered[:, 1]))
    return index


def blacklist_bin_mask(bins, blacklist_path: str | Path) -> np.ndarray:
    """True for bins overlapping at least one BED interval."""
    index = _read_blacklist(blacklist_path)
    mask = np.zeros(len(bins), dtype=bool)
    chroms = bins["chrom"].astype(str).to_numpy()
    starts = bins["start"].to_numpy(dtype=np.int64)
    ends = bins["end"].to_numpy(dtype=np.int64)
    for chrom in np.unique(chroms):
        if chrom not in index:
            continue
        rows = np.flatnonzero(chroms == chrom)
        bl_starts, prefix_max_end = index[chrom]
        positions = np.searchsorted(bl_starts, ends[rows], side="left") - 1
        valid = positions >= 0
        mask[rows[valid]] = prefix_max_end[positions[valid]] > starts[rows[valid]]
    return mask


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    source = cooler.Cooler(str(snakemake.input.cool))
    bins = source.bins()[:][["chrom", "start", "end"]]
    blocked = blacklist_bin_mask(bins, snakemake.input.blacklist)
    keep_bin = ~blocked
    nnz = int(source.info["nnz"])
    selector = source.pixels()
    chunk_size = 1_000_000

    input_contacts = 0
    retained_contacts = 0
    retained_pixels = 0

    def filtered_chunks():
        nonlocal input_contacts, retained_contacts, retained_pixels
        for lo in range(0, nnz, chunk_size):
            table = selector[lo:min(lo + chunk_size, nnz)]
            input_contacts += int(table["count"].sum())
            keep = (
                keep_bin[table["bin1_id"].to_numpy(dtype=np.int64)]
                & keep_bin[table["bin2_id"].to_numpy(dtype=np.int64)]
            )
            out = table.loc[keep, ["bin1_id", "bin2_id", "count"]].copy()
            retained_contacts += int(out["count"].sum())
            retained_pixels += len(out)
            if not out.empty:
                yield out

    destination = Path(snakemake.output.cool)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.{uuid.uuid4().hex}.tmp")
    try:
        cooler.create_cooler(
            str(temporary), bins, filtered_chunks(), ordered=True,
            symmetric_upper=True,
            metadata={
                "assembly": str(snakemake.params.assembly),
                "blacklist_sha256": _sha256(snakemake.input.blacklist),
            },
        )
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)

    write_json(
        {
            "schema": "oracle-hichip-blacklist-filter-v1",
            "sample": snakemake.wildcards.sample,
            "blacklist": str(snakemake.input.blacklist),
            "blacklist_sha256": _sha256(snakemake.input.blacklist),
            "n_bins": len(bins),
            "n_blacklisted_bins": int(blocked.sum()),
            "input_contacts": input_contacts,
            "retained_contacts": retained_contacts,
            "removed_contacts": input_contacts - retained_contacts,
            "retained_pixels": retained_pixels,
        },
        snakemake.output.json,
    )


if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
