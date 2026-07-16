"""Correctness and scale guards for consensus-loop seed de-overlap."""
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pytest


SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))
SPEC = importlib.util.spec_from_file_location(
    "consensus_performance", SCRIPTS / "build_consensus_loops.py"
)
consensus = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(consensus)


def _brute_accept(seeds, bin_size, tolerance_bins):
    radius = 2 * tolerance_bins * bin_size
    selected = []
    for seed in seeds:
        if any(
            abs(seed[0] - other[0]) <= radius
            and abs(seed[1] - other[1]) <= radius
            for other in selected
        ):
            continue
        selected.append(seed)
    return selected


def _indexed_accept(seeds, bin_size, tolerance_bins):
    index = consensus._SelectedSeedSpatialIndex(bin_size, tolerance_bins)
    selected = []
    for seed in seeds:
        if index.overlaps(seed):
            continue
        index.add(seed)
        selected.append(seed)
    return selected


@pytest.mark.parametrize("tolerance_bins", [0, 1, 2, 4])
def test_spatial_seed_index_matches_brute_force_randomized(tolerance_bins):
    rng = np.random.default_rng(20260716 + tolerance_bins)
    bin_size = 5_000
    for _ in range(40):
        grid = rng.integers(0, 80, size=(250, 2), endpoint=False)
        seeds = [tuple(map(int, row * bin_size)) for row in grid]
        assert _indexed_accept(seeds, bin_size, tolerance_bins) == _brute_accept(
            seeds, bin_size, tolerance_bins
        )


def test_spatial_seed_index_does_not_expand_conflicts_transitively():
    bin_size = 5_000
    # At tolerance=1 the exclusion radius is two bins: A conflicts with B and B
    # with C, while A does not conflict with C. B is rejected; C must still emit.
    seeds = [(0, 100_000), (10_000, 110_000), (20_000, 120_000)]
    assert _indexed_accept(seeds, bin_size, 1) == [seeds[0], seeds[2]]


def test_spatial_seed_index_scales_to_real_call_volume():
    bin_size = 5_000
    n_seeds = 75_000
    index = consensus._SelectedSeedSpatialIndex(bin_size, tolerance_bins=1)
    started = time.perf_counter()
    for i in range(n_seeds):
        # Three-bin spacing is just beyond the two-bin exclusion radius, so every
        # seed is retained and the index grows to the real upper call-set scale.
        seed = (i * 3 * bin_size, i * 3 * bin_size)
        assert not index.overlaps(seed)
        index.add(seed)
    elapsed = time.perf_counter() - started
    assert len(index) == n_seeds
    assert elapsed < 5.0, f"75k spatial lookups took {elapsed:.2f}s"
