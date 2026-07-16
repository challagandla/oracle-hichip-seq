"""
ORACLE HiChIP — shared utility functions used across scripts.
"""
import gzip
import json
import logging
from pathlib import Path
from typing import IO, TYPE_CHECKING, Iterable

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger(__name__)


GZIP_MAGIC = b"\x1f\x8b"


def open_text_auto(
    path: str | Path,
    mode: str = "rt",
    *,
    encoding: str = "utf-8",
    errors: str | None = None,
) -> IO[str]:
    """Open plain text or gzip/BGZF text by file signature, not filename.

    Public reference downloads are not consistent about retaining ``.gz`` in
    their local filename.  Conversely, users sometimes decompress a resource
    without renaming it.  Inspecting the first two bytes makes both cases safe.
    This helper is intentionally read-only so a filename cannot silently choose
    the output compression format.
    """
    if mode not in {"r", "rt"}:
        raise ValueError("open_text_auto supports read-only text mode")
    source = Path(path)
    with source.open("rb") as raw:
        compressed = raw.read(2) == GZIP_MAGIC
    kwargs: dict[str, str] = {"encoding": encoding}
    if errors is not None:
        kwargs["errors"] = errors
    if compressed:
        return gzip.open(source, "rt", **kwargs)
    return source.open("rt", **kwargs)


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
    with open_text_auto(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            chrom, size = line.strip().split("\t")[:2]
            sizes[chrom] = int(size)
    return sizes


def select_insulation_column(df: "pd.DataFrame") -> str:
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


def load_loops_bedpe(path: str | Path) -> "pd.DataFrame":
    """
    Load a FitHiChIP/generic/annotated BEDPE and normalise its core columns.

    FitHiChIP uses ``chr1/s1/e1/.../cc/P-Value_Bias/Q-Value_Bias`` rather than
    conventional BEDPE names. In particular, the raw P value and adjusted Q
    value are separate columns; preserving that distinction is essential when
    the table is later exported as graph edge attributes. Headerless 6+
    column BEDPE is also tolerated. Returns at least the six BEDPE coordinate
    columns and, when present, canonical ``score``, ``pvalue`` and ``fdr``
    columns.
    """
    # Keep the generic text/logging helpers importable in small tool-specific
    # environments (for example MACS3) that intentionally do not ship pandas.
    import pandas as pd

    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame(columns=["chrom1", "start1", "end1", "chrom2", "start2", "end2", "score", "fdr"])

    first = ""
    with open_text_auto(p, errors="ignore") as fh:
        for line in fh:
            if line.strip() and not line.lstrip().startswith("#"):
                first = line.rstrip("\n\r")
                break

    first_cols = [c.strip().lower().lstrip("#") for c in first.split("\t")]
    header_prefixes = {
        ("chrom1", "start1", "end1", "chrom2", "start2", "end2"),
        ("chr1", "s1", "e1", "chr2", "s2", "e2"),
    }
    has_header = tuple(first_cols[:6]) in header_prefixes
    if has_header:
        with open_text_auto(path) as handle:
            df = pd.read_csv(handle, sep="\t", comment="#")

        aliases = {
            "chrom1": ("chrom1", "chr1"),
            "start1": ("start1", "s1"),
            "end1": ("end1", "e1"),
            "chrom2": ("chrom2", "chr2"),
            "start2": ("start2", "s2"),
            "end2": ("end2", "e2"),
            "score": ("score", "cc"),
            "pvalue": ("pvalue", "p-value_bias", "p_value_bias", "p-value"),
            "fdr": ("fdr", "q-value_bias", "q_value_bias", "q-value"),
        }
        by_lower = {str(c).strip().lower(): c for c in df.columns}
        rename: dict[object, str] = {}
        for canonical, candidates in aliases.items():
            for candidate in candidates:
                original = by_lower.get(candidate)
                if original is not None:
                    rename[original] = canonical
                    break
        df = df.rename(columns=rename)
    else:
        with open_text_auto(path) as handle:
            df = pd.read_csv(handle, sep="\t", header=None, comment="#")
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
    for c in ("score", "pvalue", "fdr"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def write_json(obj, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, default=_default_json)


def _default_json(o):
    # NumPy is needed only when a non-native scalar reaches the JSON fallback;
    # importing it lazily keeps open_text_auto/setup_logging dependency-free.
    import numpy as np

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
