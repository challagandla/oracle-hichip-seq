"""Run Mustache only when its required balanced matrix is scientifically valid."""
import os
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from balance_utils import load_balance_report, resolution_balance  # noqa: E402
from mustache_runtime import (  # noqa: E402
    EXECUTION_BACKEND,
    MUSTACHE_HEADER,
    run_mustache_threaded,
    validate_and_sort_output,
)
from utils import load_loops_bedpe, setup_logging, write_json  # noqa: E402


LoopKey = tuple[str, int, str, int]


def canonical_loop_keys(path: str | Path, resolution: int) -> set[LoopKey]:
    """Project loop midpoints to one grid for a descriptive caller comparison."""
    if int(resolution) <= 0:
        raise ValueError("caller comparison resolution must be positive")
    loops = load_loops_bedpe(path)
    keys: set[LoopKey] = set()
    for row in loops.itertuples(index=False):
        chrom1, chrom2 = str(row.chrom1), str(row.chrom2)
        if chrom1 != chrom2:
            continue
        start1 = ((int(row.start1) + int(row.end1)) // 2 // resolution) * resolution
        start2 = ((int(row.start2) + int(row.end2)) // 2 // resolution) * resolution
        left = (chrom1, start1)
        right = (chrom2, start2)
        if left == right:
            continue
        if right < left:
            left, right = right, left
        keys.add((left[0], left[1], right[0], right[1]))
    return keys


def reciprocal_anchor_matches(
    primary: set[LoopKey],
    mustache: set[LoopKey],
    resolution: int,
    tolerance_bins: int,
) -> tuple[list[tuple[LoopKey, LoopKey]], dict[str, int]]:
    """Match reciprocal anchors with bounded grid lookups and legacy tie ranking.

    Each primary loop can only match Mustache keys in the square of
    ``(2*tolerance_bins + 1)^2`` neighboring anchor-offset combinations.  Sorting
    ``(Manhattan distance, primary key, Mustache key)`` exactly preserves the
    previous greedy one-to-one ranking without enumerating the Cartesian product.
    """
    resolution = int(resolution)
    if resolution <= 0:
        raise ValueError("caller comparison resolution must be positive")
    if tolerance_bins < 0:
        raise ValueError("caller comparison tolerance_bins must be non-negative")
    tolerance = int(tolerance_bins)
    offsets = range(-tolerance, tolerance + 1)
    edges: list[tuple[int, LoopKey, LoopKey]] = []
    grid_lookups = 0
    for primary_key in primary:
        for offset1 in offsets:
            for offset2 in offsets:
                candidate = (
                    primary_key[0],
                    primary_key[1] + offset1 * resolution,
                    primary_key[2],
                    primary_key[3] + offset2 * resolution,
                )
                grid_lookups += 1
                if candidate in mustache:
                    distance = (abs(offset1) + abs(offset2)) * resolution
                    edges.append((distance, primary_key, candidate))

    used_primary: set[LoopKey] = set()
    used_mustache: set[LoopKey] = set()
    matches: list[tuple[LoopKey, LoopKey]] = []
    for _distance, primary_key, mustache_key in sorted(edges):
        if primary_key in used_primary or mustache_key in used_mustache:
            continue
        used_primary.add(primary_key)
        used_mustache.add(mustache_key)
        matches.append((primary_key, mustache_key))
    return matches, {
        "grid_lookups": grid_lookups,
        "candidate_edges": len(edges),
    }


def caller_concordance(
    primary_path: str | Path,
    mustache_path: str | Path,
    resolution: int,
    tolerance_bins: int = 1,
) -> dict:
    """One-to-one reciprocal-anchor overlap with a documented grid tolerance."""
    if tolerance_bins < 0:
        raise ValueError("caller comparison tolerance_bins must be non-negative")
    primary = canonical_loop_keys(primary_path, resolution)
    mustache = canonical_loop_keys(mustache_path, resolution)
    matches, diagnostics = reciprocal_anchor_matches(
        primary, mustache, resolution, tolerance_bins
    )
    n_overlap = len(matches)
    union = len(primary) + len(mustache) - n_overlap
    return {
        "comparison_resolution_bp": int(resolution),
        "comparison_tolerance_bins": int(tolerance_bins),
        "comparison_matching_algorithm": "bounded_grid_greedy_v1",
        "comparison_grid_lookups": diagnostics["grid_lookups"],
        "comparison_candidate_edges": diagnostics["candidate_edges"],
        "n_primary_loops_on_grid": len(primary),
        "n_mustache_loops_on_grid": len(mustache),
        "n_overlapping_loops_on_grid": n_overlap,
        "caller_jaccard": n_overlap / union if union else None,
        "primary_supported_by_mustache_fraction": (
            n_overlap / len(primary) if primary else None
        ),
        "mustache_supported_by_primary_fraction": (
            n_overlap / len(mustache) if mustache else None
        ),
        "comparison_is_gate": False,
    }


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    resolution = int(snakemake.params.res)
    decision = resolution_balance(
        load_balance_report(snakemake.input.balance), resolution
    )
    output = Path(snakemake.output.tsv)
    output.parent.mkdir(parents=True, exist_ok=True)
    status = {
        "schema": "oracle-hichip-mustache-status-v1",
        "sample": snakemake.wildcards.sample,
        "resolution_bp": resolution,
        "balance_status": decision["status"],
        "balance_converged": decision.get("converged"),
        "execution_backend": "not_run_balance_gate",
    }
    if not decision["use_balanced"]:
        output.write_text(MUSTACHE_HEADER)
        status.update({
            "status": "NOT_ASSESSED",
            "available": False,
            "normalization": None,
            "reason": (
                "Mustache 1.3.3 requires balanced bins/weight, but ICE did not "
                "converge at the configured resolution"
            ),
            "comparison_resolution_bp": resolution,
            "comparison_tolerance_bins": int(
                snakemake.params.comparison_tolerance_bins
            ),
            "comparison_matching_algorithm": "bounded_grid_greedy_v1",
            "comparison_grid_lookups": None,
            "comparison_candidate_edges": None,
            "n_primary_loops_on_grid": len(
                canonical_loop_keys(snakemake.input.primary, resolution)
            ),
            "n_mustache_loops_on_grid": None,
            "n_overlapping_loops_on_grid": None,
            "caller_jaccard": None,
            "primary_supported_by_mustache_fraction": None,
            "mustache_supported_by_primary_fraction": None,
            "comparison_is_gate": False,
        })
        write_json(status, snakemake.output.status)
        return

    temporary_output = output.with_name(
        f".{output.name}.mustache-{os.getpid()}.tmp"
    )
    temporary_output.unlink(missing_ok=True)
    command = [
        "-f", str(snakemake.input.mcool),
        "-r", str(resolution), "-p", str(int(snakemake.threads)),
        "-o", str(temporary_output),
    ]
    try:
        with Path(snakemake.log[0]).open("a") as log:
            with redirect_stdout(log), redirect_stderr(log):
                run_mustache_threaded(command)
        n_output_rows = validate_and_sort_output(temporary_output, resolution)
        temporary_output.replace(output)
    finally:
        temporary_output.unlink(missing_ok=True)
    status.update({
        "status": "PASS",
        "available": True,
        "normalization": "ICE-balanced",
        "reason": None,
        "execution_backend": EXECUTION_BACKEND,
        "n_output_rows": n_output_rows,
    })
    status.update(
        caller_concordance(
            snakemake.input.primary,
            output,
            resolution,
            int(snakemake.params.comparison_tolerance_bins),
        )
    )
    write_json(status, snakemake.output.status)


if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
