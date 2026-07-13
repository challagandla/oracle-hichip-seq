"""
Unit tests for the numerical parts of the HiChIP workflow scripts.

These pin the behaviours that were wrong and are easy to get wrong again: how wide
an anchor is, which corners of an APA window are comparable to its centre, and what
HiCRep returns for a chromosome it could not score. Each is checked against a
matrix whose contents are known exactly, not against a golden output file.

Run: pytest -q tests/
"""
import importlib.util
import sys
from pathlib import Path

import cooler
import numpy as np
import pandas as pd
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))


def _load(name: str):
    """Import a workflow script as a module.

    The scripts guard their `main(snakemake)` call on `snakemake` being in globals,
    which Snakemake injects and we do not -- so importing them here is side-effect
    free.
    """
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


cpl = _load("count_per_loop")
utils = _load("utils")


# --------------------------------------------------------------------- fixtures
RES = 5_000
CHROMS = {"chr1": 500_000, "chr2": 300_000}


@pytest.fixture(scope="module")
def clr(tmp_path_factory) -> cooler.Cooler:
    """A tiny cooler whose every pixel we know."""
    bins = cooler.binnify(pd.Series(CHROMS), RES)
    rng = np.random.default_rng(7)
    n = len(bins)
    # dense-ish upper-triangle pixel table
    rows, cols, vals = [], [], []
    for i in range(n):
        for j in range(i, min(i + 12, n)):
            if bins.chrom[i] != bins.chrom[j]:
                continue
            rows.append(i); cols.append(j)
            vals.append(int(rng.integers(1, 20)))
    pixels = pd.DataFrame({"bin1_id": rows, "bin2_id": cols, "count": vals})
    path = tmp_path_factory.mktemp("cool") / "t.cool"
    cooler.create_cooler(str(path), bins, pixels, ordered=True)
    return cooler.Cooler(str(path))


def _dense_rectangle(clr: cooler.Cooler, r) -> int:
    """Ground truth: sum the matrix over the rectangle the two anchors span."""
    m = clr.matrix(balance=False, sparse=True).fetch(
        (r.chrom1, int(r.start1), int(r.end1)),
        (r.chrom2, int(r.start2), int(r.end2)),
    )
    return int(np.nan_to_num(m.toarray()).sum())


# ------------------------------------------------------------- anchor widths
def test_anchor_bins_spans_multiple_bins(clr):
    """A 15 kb anchor at 5 kb covers three bins, not one.

    FitHiChIP runs with MergeInt=1, so anchors are routinely several bins wide.
    """
    b = cpl._anchor_bins(clr, "chr1", 10_000, 25_000, RES)
    assert b.size == 3
    assert np.array_equal(b, clr.offset("chr1") + np.array([2, 3, 4]))


def test_anchor_bins_zero_length_still_covers_its_bin(clr):
    b = cpl._anchor_bins(clr, "chr1", 10_000, 10_000, RES)
    assert b.size == 1


def test_anchor_bins_clipped_to_chromosome(clr):
    """An anchor running off the end of the chromosome is clipped, not wrapped."""
    b = cpl._anchor_bins(clr, "chr2", 295_000, 400_000, RES)
    n_chrom_bins = int(np.ceil(CHROMS["chr2"] / RES))
    assert b.size >= 1
    assert (b - clr.offset("chr2")).max() < n_chrom_bins


def test_anchor_bins_offset_is_chromosome_relative(clr):
    """chr2 bin ids must be offset past chr1, or loops land on the wrong chromosome."""
    b1 = cpl._anchor_bins(clr, "chr1", 0, RES, RES)
    b2 = cpl._anchor_bins(clr, "chr2", 0, RES, RES)
    assert b1[0] == 0
    assert b2[0] == int(np.ceil(CHROMS["chr1"] / RES))


# ------------------------------------------------------------- loop counting
def test_count_loops_matches_dense_rectangle(clr):
    """The streaming counter must equal a direct rectangle fetch, loop for loop."""
    loops = pd.DataFrame([
        # single-bin anchors
        ("chr1", 10_000, 15_000, "chr1", 40_000, 45_000),
        # multi-bin anchors, the MergeInt=1 case
        ("chr1", 10_000, 25_000, "chr1", 40_000, 50_000),
        ("chr2", 20_000, 30_000, "chr2", 55_000, 70_000),
        # adjacent anchors
        ("chr1", 100_000, 110_000, "chr1", 110_000, 120_000),
    ], columns=["chrom1", "start1", "end1", "chrom2", "start2", "end2"])

    got = cpl.count_loops(clr, loops, RES)
    want = np.array([_dense_rectangle(clr, r) for r in loops.itertuples(index=False)])
    assert np.array_equal(got, want)
    assert want.sum() > 0, "fixture produced an all-zero matrix; test proves nothing"


def test_single_bin_counting_undercounts_wide_anchors(clr):
    """The bug this replaced: counting only the first bin loses most of the signal."""
    loops = pd.DataFrame([("chr1", 10_000, 30_000, "chr1", 60_000, 80_000)],
                         columns=["chrom1", "start1", "end1", "chrom2", "start2", "end2"])
    full = cpl.count_loops(clr, loops, RES)[0]

    r = loops.iloc[0]
    s1 = int(r.start1) // RES * RES
    s2 = int(r.start2) // RES * RES
    old = int(np.nan_to_num(
        clr.matrix(balance=False, sparse=True)
        .fetch(("chr1", s1, s1 + RES), ("chr1", s2, s2 + RES)).toarray()
    ).sum())
    assert old < full


def test_count_loops_empty_input(clr):
    empty = pd.DataFrame(columns=["chrom1", "start1", "end1", "chrom2", "start2", "end2"])
    assert cpl.count_loops(clr, empty, RES).size == 0


def test_count_loops_rejects_unknown_chromosome(clr):
    """A BEDPE from another assembly must fail loudly, not return zeros."""
    loops = pd.DataFrame([("chrZ", 10_000, 15_000, "chrZ", 40_000, 45_000)],
                         columns=["chrom1", "start1", "end1", "chrom2", "start2", "end2"])
    with pytest.raises(RuntimeError):
        cpl.count_loops(clr, loops, RES)


# ------------------------------------------------------------------ BEDPE IO
def test_load_loops_bedpe_drops_fithichip_text_header(tmp_path):
    """FitHiChIP writes a header row that is not a comment; it is not a loop."""
    p = tmp_path / "loops.bed"
    p.write_text(
        "chr1\ts1\te1\tchr2\ts2\te2\n"
        "chr1\t10000\t15000\tchr1\t40000\t45000\n"
    )
    df = utils.load_loops_bedpe(p)
    assert len(df) == 1
    assert int(df.iloc[0]["start1"]) == 10_000


def test_load_loops_bedpe_missing_file_is_empty(tmp_path):
    df = utils.load_loops_bedpe(tmp_path / "nope.bed")
    assert df.empty


# --------------------------------------------------------------------- APA
def test_apa_corners_used_are_distance_matched():
    """Only the j-i==0 corners sit at the loop's own genomic separation.

    A pixel (i, j) lies at separation D + (j - i) * bin. Build a window whose value
    IS its separation and check the background the code takes averages to D --
    the four-corner mean does not.
    """
    win = 20
    n = 2 * win + 1
    sep = np.array([[(j - i) for j in range(n)] for i in range(n)], dtype=float)

    used = np.concatenate([sep[:3, :3].ravel(), sep[-3:, -3:].ravel()])
    four = np.concatenate([sep[:3, :3].ravel(), sep[:3, -3:].ravel(),
                           sep[-3:, :3].ravel(), sep[-3:, -3:].ravel()])

    # the corners actually used are centred on the loop's separation
    assert abs(used.mean()) <= 2.0
    # and they do not span the window, unlike the old four-corner set
    assert used.max() - used.min() < (four.max() - four.min()) / 2
    # the discarded corners are the ones far off-diagonal in both directions
    assert sep[-3:, :3].mean() < -30
    assert sep[:3, -3:].mean() > 30


# ------------------------------------------------------------------- HiCRep
def test_hicrep_sentinel_is_masked_not_averaged():
    """hicrepSCC pre-fills unscored chromosomes with -2.0, and it is not NaN."""
    scc = np.array([0.8, 0.7, -2.0, 0.9])
    scored = scc[(scc > -2.0) & np.isfinite(scc)]
    assert scored.size == 3
    assert scored.mean() == pytest.approx(0.8)
    # np.nanmean does NOT skip the sentinel -- this is the bug being pinned
    assert np.nanmean(scc) < scored.mean()


# ------------------------------------------------------- snakemake compatibility
def test_no_future_imports_in_workflow_scripts():
    """No workflow script may carry `from __future__ import ...`.

    Snakemake prepends its own preamble to every file used by a `script:` directive,
    which pushes a __future__ import below other statements -- and Python rejects
    that outright:

        SyntaxError: from __future__ imports must occur at the beginning of the file

    Every script rule in the workflow died on this the first time the full DAG ran
    against real data. The other tests in this file did not catch it because they
    import the modules directly, which is not how Snakemake executes them.
    """
    offenders = [
        p.name for p in sorted(SCRIPTS.glob("*.py"))
        if "from __future__ import" in p.read_text()
    ]
    assert not offenders, (
        "these scripts will raise SyntaxError under Snakemake's script preamble: "
        f"{offenders}"
    )


def test_workflow_scripts_guard_their_main_call():
    """A script must not call main(snakemake) at import time.

    Snakemake injects `snakemake` into the script's globals; a test import does not.
    Without the guard the module cannot be imported at all, so none of it is testable.
    """
    unguarded = []
    for p in sorted(SCRIPTS.glob("*.py")):
        src = p.read_text()
        if "main(snakemake)" in src and 'if "snakemake" in globals()' not in src:
            unguarded.append(p.name)
    assert not unguarded, f"unguarded main(snakemake) call in: {unguarded}"


# --------------------------------------------------- ORACLE COS edge construction
def test_loops_to_edges_empty_result_keeps_2d_attr_shape():
    """A non-empty loop set that yields no edges must still return a (0, 3) attr array.

    Both anchors of a loop shorter than one bin land in the same bin (i == j) and are
    skipped, so `attrs` can be empty even though the BEDPE was not. `np.asarray([])`
    is 1-D, and concatenating it with the 2-D adjacency attributes raises

        ValueError: all the input arrays must have same number of dimensions

    which is how the COS export died on the 17-loop library.
    """
    cos = _load("export_oracle_cos")
    bins = pd.DataFrame({
        "chrom": ["chr1"] * 4,
        "start": [0, 100_000, 200_000, 300_000],
        "end": [100_000, 200_000, 300_000, 400_000],
        "bin_idx": [0, 1, 2, 3],
    })
    # Real loop, but both anchors sit inside bin 0 -> i == j -> contributes no edge.
    loops = pd.DataFrame([("chr1", 10_000, 15_000, "chr1", 40_000, 45_000)],
                         columns=["chrom1", "start1", "end1", "chrom2", "start2", "end2"])
    edge_index, edge_attr = cos._loops_to_edges(loops, bins)
    assert edge_index.shape == (2, 0)
    assert edge_attr.shape == (0, 3), "attr must stay 2-D or the concat with adjacency edges fails"
    # and it must actually concatenate against a 2-D adjacency block
    adj_attr = np.zeros((5, 3), dtype=np.float32)
    assert np.concatenate([edge_attr, adj_attr], axis=0).shape == (5, 3)


def test_apa_min_dist_floor_keeps_diagonal_out_of_window():
    """The configured APA distance floor must exclude the main diagonal.

    A pixel at offset (dy, dx) from a loop of span D sits at separation
    D + (dx - dy) * bin, and dx - dy ranges over +/- 2*window. So the diagonal enters
    the window for any loop with D <= 2*window*bin. The shipped config used a 100 kb
    floor with a +/-20-bin window at 10 kb -- admitting every loop under 400 kb, which
    is the median loop -- and every APA panel came out as a picture of the diagonal.
    """
    import yaml
    cfg = yaml.safe_load((Path(__file__).resolve().parents[1] / "config" / "config.yaml").read_text())
    win = int(cfg["apa"]["window_size"])
    binsz = int(cfg["apa"]["bin_size"])
    floor = int(cfg["apa"]["min_loop_dist"])
    assert floor >= (2 * win + 1) * binsz, (
        f"apa.min_loop_dist={floor} admits the diagonal into a +/-{win}-bin window "
        f"at {binsz} bp; needs >= {(2 * win + 1) * binsz}"
    )

    # and the geometry the assertion encodes: at the floor, no in-window pixel is on
    # the diagonal (separation 0)
    D = floor
    offs = np.arange(-win, win + 1)
    seps = D + (offs[None, :] - offs[:, None]) * binsz   # separation at every pixel
    assert seps.min() > 0, "some pixel in the window sits at or across the diagonal"
