"""Retain assay-stratum anchor bases supported by independent libraries."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils import setup_logging  # noqa: E402


def consensus_peaks(files: list[str], min_support: int) -> pd.DataFrame:
    if min_support < 1:
        raise ValueError("min_support must be at least one")
    records = []
    for path in files:
        sample = Path(path).name.removesuffix("_peaks.bed")
        peaks = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2],
                            names=["chrom", "start", "end"])
        peaks[["start", "end"]] = peaks[["start", "end"]].apply(
            pd.to_numeric, errors="raise"
        ).astype(int)
        if (peaks["end"] <= peaks["start"]).any():
            raise ValueError(f"peak file contains end <= start intervals: {path}")
        peaks["sample"] = sample
        records.extend(peaks.itertuples(index=False, name=None))
    if not records:
        return pd.DataFrame(
            columns=["chrom", "start", "end", "sample_support", "support_samples"]
        )

    supported = []
    table = pd.DataFrame(
        records, columns=["chrom", "start", "end", "sample"]
    )
    for chrom, chrom_group in table.groupby("chrom", sort=False):
        # First merge within each sample. A sample with overlapping peak records
        # still contributes one unit of biological support, never two.
        merged_by_sample: list[tuple[int, int, str]] = []
        for sample, group in chrom_group.groupby("sample", sort=False):
            current_start = current_end = None
            for row in group.sort_values(["start", "end"], kind="stable").itertuples(
                index=False
            ):
                start, end = int(row.start), int(row.end)
                if current_end is None or start >= current_end:
                    if current_end is not None:
                        merged_by_sample.append(
                            (int(current_start), int(current_end), str(sample))
                        )
                    current_start, current_end = start, end
                else:
                    current_end = max(int(current_end), end)
            if current_end is not None:
                merged_by_sample.append(
                    (int(current_start), int(current_end), str(sample))
                )

        events: dict[int, dict[str, set[str]]] = {}
        for start, end, sample in merged_by_sample:
            events.setdefault(start, {"add": set(), "remove": set()})["add"].add(sample)
            events.setdefault(end, {"add": set(), "remove": set()})["remove"].add(sample)

        active: set[str] = set()
        previous = None
        for position in sorted(events):
            if previous is not None and position > previous and len(active) >= min_support:
                supported.append(
                    (chrom, previous, position, len(active), ",".join(sorted(active)))
                )
            active.difference_update(events[position]["remove"])
            active.update(events[position]["add"])
            previous = position

    return pd.DataFrame(
        supported,
        columns=["chrom", "start", "end", "sample_support", "support_samples"],
    )


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    consensus = consensus_peaks(list(snakemake.input.peaks), int(snakemake.params.min_support))
    if consensus.empty and not bool(getattr(snakemake.params, "allow_empty", False)):
        raise RuntimeError("No independently supported peaks remain for this assay stratum")
    Path(snakemake.output.bed).parent.mkdir(parents=True, exist_ok=True)
    consensus[["chrom", "start", "end"]].to_csv(
        snakemake.output.bed, sep="\t", index=False, header=False
    )
    consensus.to_csv(snakemake.output.audit, sep="\t", index=False)


if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
