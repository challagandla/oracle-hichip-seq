"""Small, dependency-free helpers for configured locus viewpoints."""
from typing import Any, Mapping


def resolve_viewpoint(region: Mapping[str, Any]) -> tuple[int, str]:
    """Return and validate a region's 0-based viewpoint and display label."""
    if "viewpoint" not in region:
        raise ValueError(
            f"visualisation region {region.get('name', '<unnamed>')!r} has no explicit "
            "viewpoint; add a 0-based genomic coordinate to config.yaml"
        )

    start = int(region["start"])
    end = int(region["end"])
    viewpoint = int(region["viewpoint"])
    if end <= start:
        raise ValueError(f"invalid visualisation interval: start={start}, end={end}")
    if not start <= viewpoint < end:
        raise ValueError(
            f"viewpoint {viewpoint} lies outside the configured half-open interval "
            f"[{start}, {end})"
        )

    label = str(region.get("viewpoint_label") or f"viewpoint {viewpoint:,}")
    return viewpoint, label
