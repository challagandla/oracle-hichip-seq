"""Build a fixed-bin, multi-library-supported candidate set for differential testing.

Independent FitHiChIP call sets often place the same contact one bin apart. This
module maps both anchors to the configured grid and uses a deterministic,
representative-centred reciprocal-anchor tolerance. Each library contributes at
most one vote to a candidate and a neighbourhood never expands transitively.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_loops_bedpe, open_text_auto, read_chromsizes, setup_logging  # noqa: E402

COORDS = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]


def _primary_chromosome(chrom: pd.Series) -> pd.Series:
    return chrom.astype(str).str.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2]|X)")


def canonicalise_loops(loops: pd.DataFrame, bin_size: int) -> pd.DataFrame:
    """Map each anchor midpoint to one fixed-width bin and deduplicate."""
    if loops.empty:
        return pd.DataFrame(columns=COORDS)
    out = loops[COORDS].copy()
    mid1 = (out["start1"].astype(np.int64) + out["end1"].astype(np.int64)) // 2
    mid2 = (out["start2"].astype(np.int64) + out["end2"].astype(np.int64)) // 2
    out["start1"] = (mid1 // bin_size) * bin_size
    out["end1"] = out["start1"] + bin_size
    out["start2"] = (mid2 // bin_size) * bin_size
    out["end2"] = out["start2"] + bin_size
    out = out[
        _primary_chromosome(out["chrom1"])
        & _primary_chromosome(out["chrom2"])
        & (out["chrom1"] == out["chrom2"])
    ]

    # Enforce an upper-triangle representation before dropping duplicates.
    swap = (out["chrom1"] > out["chrom2"]) | (
        (out["chrom1"] == out["chrom2"]) & (out["start1"] > out["start2"])
    )
    if swap.any():
        left = out.loc[swap, ["chrom1", "start1", "end1"]].to_numpy(copy=True)
        out.loc[swap, ["chrom1", "start1", "end1"]] = out.loc[
            swap, ["chrom2", "start2", "end2"]
        ].to_numpy()
        out.loc[swap, ["chrom2", "start2", "end2"]] = left
    out = out[out["start1"] != out["start2"]]
    return out.drop_duplicates(COORDS).reset_index(drop=True)


def _read_blacklist(path: str | Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    with open_text_auto(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end = line.rstrip().split("\t")[:3]
            by_chrom.setdefault(chrom, []).append((int(start), int(end)))
    index = {}
    for chrom, intervals in by_chrom.items():
        ordered = np.asarray(sorted(intervals), dtype=np.int64)
        index[chrom] = (ordered[:, 0], np.maximum.accumulate(ordered[:, 1]))
    return index


def _overlaps_blacklist(
    chroms: pd.Series, starts: pd.Series, ends: pd.Series,
    index: dict[str, tuple[np.ndarray, np.ndarray]],
) -> np.ndarray:
    hit = np.zeros(len(chroms), dtype=bool)
    chrom_values = chroms.astype(str).to_numpy()
    start_values = starts.to_numpy(dtype=np.int64)
    end_values = ends.to_numpy(dtype=np.int64)
    for chrom in np.unique(chrom_values):
        if chrom not in index:
            continue
        rows = np.flatnonzero(chrom_values == chrom)
        bl_starts, prefix_max_end = index[chrom]
        pos = np.searchsorted(bl_starts, end_values[rows], side="left") - 1
        valid = pos >= 0
        hit[rows[valid]] = prefix_max_end[pos[valid]] > start_values[rows[valid]]
    return hit


def _interval_overlaps(
    index: dict[str, tuple[np.ndarray, np.ndarray]],
    chrom: str,
    start: int,
    end: int,
) -> bool:
    if chrom not in index:
        return False
    starts, prefix_max_end = index[chrom]
    pos = int(np.searchsorted(starts, end, side="left") - 1)
    return pos >= 0 and int(prefix_max_end[pos]) > int(start)


class _SelectedSeedSpatialIndex:
    """Bounded lookup for already-emitted representative loop pixels.

    Candidate seeds live on a fixed two-dimensional bin grid.  Two emitted
    candidates conflict when both anchor coordinates are within
    ``2 * tolerance_bins`` grid steps.  Buckets are one grid step wider than
    that exclusion radius, so a possible conflict can only occur in the seed's
    own bucket or one of its eight neighbours.  At most one accepted seed can
    occupy a bucket, keeping lookup O(1) in the number of loop calls.
    """

    def __init__(self, bin_size: int, tolerance_bins: int) -> None:
        if int(bin_size) <= 0:
            raise ValueError("bin_size must be positive")
        if int(tolerance_bins) < 0:
            raise ValueError("tolerance_bins must be non-negative")
        self.bin_size = int(bin_size)
        self.radius = 2 * int(tolerance_bins)
        self.bucket_width = self.radius + 1
        self._buckets: dict[tuple[int, int], tuple[int, int]] = {}

    def _grid(self, seed: tuple[int, int]) -> tuple[int, int]:
        if seed[0] % self.bin_size or seed[1] % self.bin_size:
            raise ValueError(f"seed is not aligned to the {self.bin_size}-bp grid: {seed}")
        return seed[0] // self.bin_size, seed[1] // self.bin_size

    def _bucket(self, grid: tuple[int, int]) -> tuple[int, int]:
        return grid[0] // self.bucket_width, grid[1] // self.bucket_width

    def overlaps(self, seed: tuple[int, int]) -> bool:
        grid = self._grid(seed)
        bucket = self._bucket(grid)
        for delta1 in (-1, 0, 1):
            for delta2 in (-1, 0, 1):
                other = self._buckets.get((bucket[0] + delta1, bucket[1] + delta2))
                if other is not None and (
                    abs(grid[0] - other[0]) <= self.radius
                    and abs(grid[1] - other[1]) <= self.radius
                ):
                    return True
        return False

    def add(self, seed: tuple[int, int]) -> None:
        grid = self._grid(seed)
        bucket = self._bucket(grid)
        if bucket in self._buckets:
            raise ValueError(f"overlapping selected seed inserted into bucket {bucket}")
        self._buckets[bucket] = grid

    def __len__(self) -> int:
        return len(self._buckets)


def build_consensus(
    loop_files: list[str], bin_size: int, min_sample_support: int,
    blacklist: str | Path | None = None,
    required_anchor_bed: str | Path | None = None,
    tolerance_bins: int = 1,
    chromsizes: dict[str, int] | None = None,
) -> pd.DataFrame:
    observations = []
    for path in loop_files:
        canonical = canonicalise_loops(load_loops_bedpe(path), bin_size)
        canonical["sample"] = Path(path).parent.name
        observations.append(canonical)
    if not observations:
        return pd.DataFrame(columns=COORDS + ["sample_support"])
    all_calls = pd.concat(observations, ignore_index=True)
    if all_calls.empty:
        return pd.DataFrame(columns=COORDS + ["sample_support"])

    blacklist_index: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    if blacklist and Path(blacklist).exists():
        blacklist_index = _read_blacklist(blacklist)
        blocked = _overlaps_blacklist(
            all_calls["chrom1"], all_calls["start1"], all_calls["end1"], blacklist_index
        ) | _overlaps_blacklist(
            all_calls["chrom2"], all_calls["start2"], all_calls["end2"], blacklist_index
        )
        all_calls = all_calls.loc[~blocked].copy()

    # For contact-map-held-out APA, require at least one Peak-to-ALL anchor to be
    # supported by sibling donors without using the scored sample's peak calls.
    # The primary sibling loop calls still come from the cohort-wide FitHiChIP
    # search space; this filter removes target-only anchor leakage and that residual
    # search-space dependence is recorded in the APA metadata.
    if required_anchor_bed is not None:
        anchor_path = Path(required_anchor_bed)
        if not anchor_path.exists():
            raise FileNotFoundError(f"required anchor BED is missing: {anchor_path}")
        if anchor_path.stat().st_size == 0:
            all_calls = all_calls.iloc[0:0].copy()
        else:
            anchor_index = _read_blacklist(anchor_path)
            supported_anchor = _overlaps_blacklist(
                all_calls["chrom1"], all_calls["start1"], all_calls["end1"],
                anchor_index,
            ) | _overlaps_blacklist(
                all_calls["chrom2"], all_calls["start2"], all_calls["end2"],
                anchor_index,
            )
            all_calls = all_calls.loc[supported_anchor].copy()

    if tolerance_bins < 0:
        raise ValueError("candidate tolerance_bins must be non-negative")

    records: list[dict] = []
    # Calls are cis and upper-triangular after canonicalisation, so a 2-D grid of
    # anchor starts is sufficient. Candidate neighbourhoods are centred on one
    # observed call and never expanded from newly included neighbours; A~B and B~C
    # therefore cannot pull in C unless C is also within tolerance of the chosen
    # representative. The highest-support representative wins deterministically.
    for chrom, group in all_calls.groupby("chrom1", sort=True):
        observations: dict[tuple[int, int], set[str]] = {}
        for row in group.drop_duplicates(COORDS + ["sample"]).itertuples(index=False):
            observations.setdefault((int(row.start1), int(row.start2)), set()).add(
                str(row.sample)
            )

        def neighbours(seed: tuple[int, int]) -> list[tuple[int, int]]:
            return [
                coord
                for delta1 in range(-int(tolerance_bins), int(tolerance_bins) + 1)
                for delta2 in range(-int(tolerance_bins), int(tolerance_bins) + 1)
                if (
                    coord := (
                        seed[0] + delta1 * int(bin_size),
                        seed[1] + delta2 * int(bin_size),
                    )
                ) in observations
            ]

        ranked = []
        for seed in observations:
            nearby = neighbours(seed)
            samples = set().union(*(observations[coord] for coord in nearby))
            ranked.append((-len(samples), seed[0], seed[1], seed, nearby))
        ranked.sort(key=lambda item: item[:3])

        assigned: set[tuple[tuple[int, int], str]] = set()
        selected_seeds = _SelectedSeedSpatialIndex(bin_size, tolerance_bins)
        for _neg_support, _start1, _start2, seed, nearby in ranked:
            if not any(
                (seed, sample) not in assigned for sample in observations[seed]
            ):
                continue
            # Each emitted candidate is counted over a full tolerance footprint.
            # Do not emit another candidate whose two anchor footprints overlap,
            # or the same contact pixels would enter the differential matrix twice.
            if selected_seeds.overlaps(seed):
                continue
            available = {
                sample
                for coord in nearby
                for sample in observations[coord]
                if (coord, sample) not in assigned
            }
            if len(available) < min_sample_support:
                continue
            radius = int(tolerance_bins) * int(bin_size)
            chrom_end = (
                int(chromsizes[str(chrom)])
                if chromsizes is not None and str(chrom) in chromsizes
                else None
            )
            footprint = (
                max(0, seed[0] - radius),
                min(seed[0] + bin_size + radius, chrom_end)
                if chrom_end is not None else seed[0] + bin_size + radius,
                max(0, seed[1] - radius),
                min(seed[1] + bin_size + radius, chrom_end)
                if chrom_end is not None else seed[1] + bin_size + radius,
            )
            if _interval_overlaps(
                blacklist_index, str(chrom), footprint[0], footprint[1]
            ) or _interval_overlaps(
                blacklist_index, str(chrom), footprint[2], footprint[3]
            ):
                continue
            records.append({
                "chrom1": chrom,
                "start1": footprint[0],
                "end1": footprint[1],
                "chrom2": chrom,
                "start2": footprint[2],
                "end2": footprint[3],
                "representative_start1": seed[0],
                "representative_start2": seed[1],
                "sample_support": len(available),
                "support_samples": ",".join(sorted(available)),
                "tolerance_bins": int(tolerance_bins),
            })
            # Consume the representative-centred neighbourhood. This prevents a
            # sample from voting again through an overlapping shifted candidate.
            assigned.update(
                (coord, sample)
                for coord in nearby
                for sample in observations[coord]
            )
            selected_seeds.add(seed)

    columns = COORDS + [
        "representative_start1", "representative_start2", "sample_support",
        "support_samples", "tolerance_bins",
    ]
    support = pd.DataFrame(records, columns=columns)
    return support.sort_values(COORDS, kind="stable").reset_index(drop=True)


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    chromsizes = read_chromsizes(snakemake.input.chromsizes)
    consensus = build_consensus(
        list(snakemake.input.loops),
        int(snakemake.params.bin_size),
        int(snakemake.params.min_sample_support),
        snakemake.input.blacklist,
        getattr(snakemake.input, "anchors", None),
        int(getattr(snakemake.params, "tolerance_bins", 1)),
        chromsizes,
    )
    Path(snakemake.output.audit).parent.mkdir(parents=True, exist_ok=True)
    consensus.to_csv(snakemake.output.audit, sep="\t", index=False)
    consensus[COORDS].to_csv(snakemake.output.bedpe, sep="\t", index=False, header=False)
    if consensus.empty and not bool(getattr(snakemake.params, "allow_empty", False)):
        raise RuntimeError(
            "No multi-library-supported differential candidates remained; inspect the "
            "per-sample call sets or lower differential.min_sample_support explicitly."
        )


if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
