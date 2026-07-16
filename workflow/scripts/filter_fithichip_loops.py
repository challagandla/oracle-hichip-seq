"""Strictly validate and filter FitHiChIP's q-thresholded reporting calls."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from build_consensus_loops import _overlaps_blacklist, _read_blacklist
from utils import open_text_auto


RAW_COORDS = ["chr1", "s1", "e1", "chr2", "s2", "e2"]
RAW_REQUIRED = RAW_COORDS + ["cc", "P-Value_Bias", "Q-Value_Bias"]
CANONICAL_COORDS = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
PRIMARY = r"chr(?:[1-9]|1[0-9]|2[0-2]|X)"


def _header(path: str | Path) -> list[str]:
    with open_text_auto(path, errors="strict") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                return line.rstrip("\r\n").split("\t")
    raise ValueError("FitHiChIP reporting call file is empty and has no header")


def _strict_native_calls(
    source: str | Path, min_reads: int, q_threshold: float
) -> pd.DataFrame:
    if min_reads < 0:
        raise ValueError("min_reads must be non-negative")
    if not np.isfinite(q_threshold) or not 0 < q_threshold <= 1:
        raise ValueError("q_threshold must be finite and in (0,1]")
    header = _header(source)
    if "cc" not in header and min_reads > 0:
        raise ValueError(
            "FitHiChIP output has no cc/score column, so fithichip.min_reads "
            "cannot be applied"
        )
    if len(header) != len(set(header)):
        raise ValueError("FitHiChIP reporting header contains duplicate columns")
    missing = [column for column in RAW_REQUIRED if column not in header]
    if missing:
        raise ValueError(f"FitHiChIP reporting calls lack required columns: {missing}")
    if header[:6] != RAW_COORDS:
        raise ValueError(
            "FitHiChIP reporting calls must start with chr1/s1/e1/chr2/s2/e2"
        )

    with open_text_auto(source, errors="strict") as handle:
        frame = pd.read_csv(
            handle, sep="\t", comment="#", dtype=str, keep_default_na=False
        )
    numeric_columns = ["s1", "e1", "s2", "e2", "cc"]
    for column in numeric_columns:
        numeric = pd.to_numeric(frame[column], errors="coerce")
        if numeric.isna().any() or not np.equal(numeric, np.floor(numeric)).all():
            raise ValueError(f"FitHiChIP reporting column {column!r} must be integer")
        frame[column] = numeric.astype(np.int64)
    for column in ["P-Value_Bias", "Q-Value_Bias"]:
        numeric = pd.to_numeric(frame[column], errors="coerce")
        invalid = numeric.isna() | (~np.isfinite(numeric)) | (numeric < 0) | (numeric > 1)
        if invalid.any():
            examples = frame.loc[invalid, column].drop_duplicates().head(5).tolist()
            raise ValueError(
                f"FitHiChIP reporting column {column!r} must contain finite "
                f"probabilities in [0,1]; invalid values: {examples}"
            )
        frame[column] = numeric.astype(float)

    invalid_geometry = (
        (frame["s1"] < 0)
        | (frame["s2"] < 0)
        | (frame["e1"] <= frame["s1"])
        | (frame["e2"] <= frame["s2"])
        | (frame["cc"] <= 0)
    )
    if invalid_geometry.any():
        raise ValueError(
            f"FitHiChIP reporting calls contain {int(invalid_geometry.sum())} "
            "invalid coordinate/count rows"
        )
    tolerance = max(1e-12, abs(q_threshold) * 1e-9)
    above = frame["Q-Value_Bias"] > q_threshold + tolerance
    if above.any():
        raise ValueError(
            f"FitHiChIP reporting calls contain {int(above.sum())} rows above "
            f"the configured q threshold {q_threshold}"
        )
    return frame.rename(columns={
        "chr1": "chrom1", "s1": "start1", "e1": "end1",
        "chr2": "chrom2", "s2": "start2", "e2": "end2",
        "cc": "score", "P-Value_Bias": "pvalue", "Q-Value_Bias": "fdr",
    })


def filter_loops(
    source: str | Path,
    blacklist: str | Path,
    min_reads: int = 0,
    q_threshold: float = 1.0,
) -> tuple[pd.DataFrame, dict[str, int | float]]:
    loops = _strict_native_calls(source, min_reads, q_threshold)
    before = len(loops)
    primary = loops["chrom1"].astype(str).str.fullmatch(PRIMARY)
    primary &= loops["chrom2"].astype(str).str.fullmatch(PRIMARY)
    primary &= loops["chrom1"].astype(str) == loops["chrom2"].astype(str)
    loops = loops.loc[primary].copy()
    after_primary = len(loops)

    index = _read_blacklist(blacklist)
    blocked = _overlaps_blacklist(
        loops["chrom1"], loops["start1"], loops["end1"], index
    ) | _overlaps_blacklist(
        loops["chrom2"], loops["start2"], loops["end2"], index
    )
    loops = loops.loc[~blocked].copy()
    after_blacklist = len(loops)
    if min_reads > 0:
        loops = loops.loc[loops["score"] >= min_reads].copy()
    return loops, {
        "configured_q_threshold": q_threshold,
        "input": before,
        "removed_non_primary_or_trans": before - after_primary,
        "removed_blacklist": after_primary - after_blacklist,
        "removed_below_min_reads": after_blacklist - len(loops),
        "retained": len(loops),
    }


def cli() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--blacklist", required=True)
    parser.add_argument("--min-reads", type=int, default=0)
    parser.add_argument("--q-threshold", type=float, required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--audit", required=True)
    args = parser.parse_args()
    loops, audit = filter_loops(
        args.input, args.blacklist, args.min_reads, args.q_threshold
    )
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    loops.to_csv(args.output, sep="\t", index=False)
    Path(args.audit).write_text(
        "\n".join(f"{key}\t{value}" for key, value in audit.items()) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    cli()
