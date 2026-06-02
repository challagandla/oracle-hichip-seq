"""
ORACLE HiChIP — shared utility functions used across scripts.
"""
from __future__ import annotations

import gzip
import json
import logging
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def setup_logging(log_path: str | Path | None, level: int = logging.INFO) -> None:
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if log_path:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_path))
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=handlers,
        force=True,
    )


def read_chromsizes(path: str | Path) -> dict[str, int]:
    """chrom\tsize → dict."""
    sizes: dict[str, int] = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:  # type: ignore[arg-type]
        for line in fh:
            if not line.strip():
                continue
            chrom, size = line.strip().split("\t")[:2]
            sizes[chrom] = int(size)
    return sizes


def select_insulation_column(df: pd.DataFrame) -> str:
    """Return the score column from a cooltools insulation table."""
    exact_names = ("log2_insulation_score", "insulation_score", "insulation")
    for name in exact_names:
        if name in df.columns:
            return name

    excluded_tokens = ("boundary", "is_boundary", "bad_bin", "valid_pixels", "n_valid")
    prefixes = ("log2_insulation_score_", "insulation_score_", "insulation_")
    for prefix in prefixes:
        matches = [
            col for col in df.columns
            if str(col).lower().startswith(prefix)
            and not any(token in str(col).lower() for token in excluded_tokens)
        ]
        if matches:
            return matches[0]

    matches = [
        col for col in df.columns
        if "insulation" in str(col).lower()
        and not any(token in str(col).lower() for token in excluded_tokens)
    ]
    if matches:
        return matches[0]

    raise ValueError(f"Could not identify insulation score column. Columns: {list(df.columns)}")


def load_loops_bedpe(path: str | Path) -> pd.DataFrame:
    """
    Load a FitHiChIP/generic/annotated BEDPE. Headered annotated BEDPE and
    headerless 6+ column BEDPE are both tolerated. Returns at least:
    [chrom1,start1,end1,chrom2,start2,end2,score,fdr] when available.
    """
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame(columns=["chrom1", "start1", "end1", "chrom2", "start2", "end2", "score", "fdr"])

    first = p.read_text(errors="ignore").splitlines()[0]
    has_header = first.lower().startswith("chrom1\t") or first.lower().startswith("chrom1,")
    if has_header:
        df = pd.read_csv(path, sep="\t", comment="#")
    else:
        df = pd.read_csv(path, sep="\t", header=None, comment="#")
        base_cols = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "score", "fdr"]
        df = df.rename(columns={i: c for i, c in enumerate(base_cols[: df.shape[1]])})

    required = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"BEDPE {path} is missing required columns: {missing}")
    for c in ("start1", "end1", "start2", "end2"):
        df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
    df = df.dropna(subset=["start1", "end1", "start2", "end2"]).copy()
    for c in ("start1", "end1", "start2", "end2"):
        df[c] = df[c].astype(int)
    return df


def write_json(obj, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, default=_default_json)


def _default_json(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"Cannot JSON-serialise {type(o).__name__}")


def passing(value: float, *, ge: float | None = None, le: float | None = None) -> bool:
    if ge is not None and value < ge:
        return False
    if le is not None and value > le:
        return False
    return True


def chunks(it: Iterable, n: int):
    """Yield n-sized chunks from iterable."""
    buf: list = []
    for x in it:
        buf.append(x)
        if len(buf) >= n:
            yield buf
            buf = []
    if buf:
        yield buf
