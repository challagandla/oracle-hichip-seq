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
            chrom, size = line.strip().split("\t")[:2]
            sizes[chrom] = int(size)
    return sizes


def load_loops_bedpe(path: str | Path) -> pd.DataFrame:
    """
    Load a FitHiChIP / generic BEDPE. We tolerate either 6-col (BED6 BEDPE)
    or full FitHiChIP output. Returns canonical columns:
    [chrom1,start1,end1,chrom2,start2,end2,score,fdr]
    """
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    df = df.rename(columns={i: c for i, c in enumerate(
        ["chrom1", "start1", "end1", "chrom2", "start2", "end2",
         "score", "fdr"][: df.shape[1]]
    )})
    for c in ("start1", "end1", "start2", "end2"):
        if c in df.columns:
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
