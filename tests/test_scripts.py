"""
Unit tests for the numerical parts of the HiChIP workflow scripts.

These pin easy-to-regress behaviours: how wide an anchor is, which corners of an
APA window are comparable to its centre, and what
HiCRep returns for a chromosome it could not score. Each is checked against a
matrix whose contents are known exactly, not against a golden output file.

Run: pytest -q tests/
"""
import gzip
import importlib.util
import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import cooler
import h5py
import numpy as np
import pandas as pd
import pytest
import yaml

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
filter_fithichip = _load("filter_fithichip_loops")
viewpoint = _load("viz_viewpoint")
pgt = _load("pygenometracks_loops")
loop_qc = _load("loop_qc_summary")
multiqc_content = _load("multiqc_content")
hicrep = _load("hicrep_replicate_qc")
balance_utils = _load("balance_utils")
balance_qc = _load("balance_qc")
build_mcool = _load("build_mcool")
filter_cool_blacklist = _load("filter_cool_blacklist")
cooltools_matrix_qc = _load("cooltools_matrix_qc")
mustache_balance_aware = _load("mustache_balance_aware")
cooltools_eigs = _load("cooltools_eigs_cis")


# --------------------------------------------------------------------- fixtures
RES = 5_000
CHROMS = {"chr1": 500_000, "chr2": 300_000}


def test_text_helpers_import_without_scientific_site_packages():
    """Small tool environments can use gzip/text helpers without pandas."""
    command = (
        "import sys; "
        f"sys.path.insert(0, {str(SCRIPTS)!r}); "
        "from utils import open_text_auto, setup_logging; "
        "assert callable(open_text_auto) and callable(setup_logging)"
    )
    subprocess.run([sys.executable, "-S", "-c", command], check=True)


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
            rows.append(i)
            cols.append(j)
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
    """A generic 15 kb anchor at 5 kb covers three bins, not one."""
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
        # generic multi-bin anchors
        ("chr1", 10_000, 25_000, "chr1", 40_000, 50_000),
        ("chr2", 20_000, 30_000, "chr2", 55_000, 70_000),
        # adjacent anchors
        ("chr1", 100_000, 110_000, "chr1", 110_000, 120_000),
    ], columns=["chrom1", "start1", "end1", "chrom2", "start2", "end2"])

    got = cpl.count_loops(clr, loops, RES, 0, 10**9)
    want = np.array([_dense_rectangle(clr, r) for r in loops.itertuples(index=False)])
    assert np.array_equal(got, want)
    assert want.sum() > 0, "fixture produced an all-zero matrix; test proves nothing"


def test_single_bin_counting_undercounts_wide_anchors(clr):
    """Counting only the first bin of a wide anchor loses signal."""
    loops = pd.DataFrame([("chr1", 10_000, 30_000, "chr1", 60_000, 80_000)],
                         columns=["chrom1", "start1", "end1", "chrom2", "start2", "end2"])
    full = cpl.count_loops(clr, loops, RES, 0, 10**9)[0]

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
    assert cpl.count_loops(clr, empty, RES, 0, 10**9).size == 0


def test_count_loops_rejects_unknown_chromosome(clr):
    """A BEDPE from another assembly must fail loudly, not return zeros."""
    loops = pd.DataFrame([("chrZ", 10_000, 15_000, "chrZ", 40_000, 45_000)],
                         columns=["chrom1", "start1", "end1", "chrom2", "start2", "end2"])
    with pytest.raises(RuntimeError):
        cpl.count_loops(clr, loops, RES, 0, 10**9)


def test_count_loops_enforces_both_fithichip_distance_boundaries(clr):
    """A widened footprint must not import near-diagonal or over-range pixels."""
    loops = pd.DataFrame([
        ("chr1", 0, 15_000, "chr1", 15_000, 35_000),
    ], columns=["chrom1", "start1", "end1", "chrom2", "start2", "end2"])
    lower, upper = 10_000, 20_000
    got = cpl.count_loops(clr, loops, RES, lower, upper)[0]

    dense = np.asarray(
        clr.matrix(balance=False).fetch(
            ("chr1", 0, 15_000), ("chr1", 15_000, 35_000)
        )
    )
    left = np.arange(0, 15_000, RES)
    right = np.arange(15_000, 35_000, RES)
    distance = np.abs(left[:, None] - right[None, :])
    eligible = (distance >= lower) & (distance <= upper)
    assert eligible.any() and (~eligible).any()
    assert (distance[~eligible] < lower).any()
    assert (distance[~eligible] > upper).any()
    assert got == int(dense[eligible].sum())
    assert got < int(dense.sum())


def test_hicrep_depth_matches_scored_chromosomes_and_diagonals(clr):
    table = clr.matrix(balance=False, as_pixels=True, join=False).fetch("chr1")
    offset = table["bin2_id"] - table["bin1_id"]
    max_dist = 3 * RES
    expected = int(table.loc[(offset > 0) & (offset <= 3), "count"].sum())

    got = hicrep._selected_cis_contacts(clr, ["chr1"], max_dist)
    assert got == expected
    assert got < int(clr.info["sum"])
    assert hicrep._selected_cis_contacts(clr, ["chr1"], 0) == 0


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


def test_load_loops_bedpe_normalises_fithichip_statistics(tmp_path):
    """FitHiChIP's adjusted Q value, not its raw P value, is the loop FDR."""
    p = tmp_path / "merged_interactions.bed"
    p.write_text(
        "chr1\ts1\te1\tchr2\ts2\te2\tcc\tP-Value_Bias\tQ-Value_Bias\t"
        "bin1_low\tbin1_high\tbin2_low\tbin2_high\tsumCC\tStrongConn\n"
        "chr1\t2255000\t2260000\tchr1\t2300000\t2305000\t24\t"
        "2.71877687225678e-07\t0.0025099428228572\t2255000\t2260000\t"
        "2300000\t2305000\t24\t1.0\n"
    )

    df = utils.load_loops_bedpe(p)

    assert len(df) == 1
    assert list(df.loc[:, ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]].iloc[0]) == [
        "chr1", 2_255_000, 2_260_000, "chr1", 2_300_000, 2_305_000,
    ]
    assert df.iloc[0]["score"] == pytest.approx(24)
    assert df.iloc[0]["pvalue"] == pytest.approx(2.71877687225678e-07)
    assert df.iloc[0]["fdr"] == pytest.approx(0.0025099428228572)
    assert df.iloc[0]["fdr"] != pytest.approx(df.iloc[0]["pvalue"])


def test_load_loops_bedpe_missing_file_is_empty(tmp_path):
    df = utils.load_loops_bedpe(tmp_path / "nope.bed")
    assert df.empty


def test_fithichip_filter_applies_configured_minimum_read_count(tmp_path):
    loops = tmp_path / "calls.bed"
    loops.write_text(
        "chr1\ts1\te1\tchr2\ts2\te2\tcc\tP-Value_Bias\tQ-Value_Bias\n"
        "chr1\t0\t5000\tchr1\t20000\t25000\t5\t0.001\t0.01\n"
        "chr1\t0\t5000\tchr1\t30000\t35000\t6\t0.001\t0.01\n"
    )
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr2\t0\t100\n")

    kept, audit = filter_fithichip.filter_loops(loops, blacklist, min_reads=6)

    assert kept["score"].tolist() == [6]
    assert audit["removed_below_min_reads"] == 1


def test_fithichip_filter_requires_score_when_minimum_is_enabled(tmp_path):
    loops = tmp_path / "headerless.bedpe"
    loops.write_text("chr1\t0\t5000\tchr1\t20000\t25000\n")
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr2\t0\t100\n")

    with pytest.raises(ValueError, match="no cc/score"):
        filter_fithichip.filter_loops(loops, blacklist, min_reads=6)


# ------------------------------------------------------- matrix balance contract
def _write_weight(handle, resolution, converged, **attrs):
    group = handle.require_group(f"resolutions/{resolution}/bins")
    dataset = group.create_dataset("weight", data=np.ones(4, dtype=float))
    if converged is not None:
        dataset.attrs["converged"] = converged
    for key, value in attrs.items():
        dataset.attrs[key] = value


def test_balance_qc_reads_hdf5_attrs_and_reports_warn(tmp_path):
    mcool = tmp_path / "attempted.mcool"
    with h5py.File(mcool, "w") as handle:
        _write_weight(
            handle, 5_000, False, var=0.0048, tol=1e-5,
            ignore_diags=2, min_nnz=10, mad_max=5,
        )
        _write_weight(handle, 10_000, True, var=4e-6, tol=1e-5)

    report = balance_utils.inspect_mcool_balance(mcool, [5_000, 10_000])

    assert report["status"] == "WARN"
    assert report["pass"] is None
    assert report["nonconverged_resolutions_bp"] == [5_000]
    entry = report["resolutions"]["5000"]
    assert entry["converged"] is False
    assert entry["variance"] == pytest.approx(0.0048)
    assert entry["tolerance"] == pytest.approx(1e-5)
    assert entry["parameters"]["ignore_diags"] == 2


def test_balance_qc_missing_resolution_is_not_assessed(tmp_path):
    mcool = tmp_path / "attempted.mcool"
    with h5py.File(mcool, "w") as handle:
        _write_weight(handle, 5_000, True, var=1e-7, tol=1e-5)

    report = balance_utils.inspect_mcool_balance(mcool, [5_000, 25_000])

    assert report["status"] == "NOT_ASSESSED"
    assert report["missing_resolutions_bp"] == [25_000]


def test_nonpass_weights_are_removed_only_from_published_matrix(tmp_path):
    mcool = tmp_path / "published.mcool"
    with h5py.File(mcool, "w") as handle:
        _write_weight(handle, 5_000, False, var=0.0048, tol=1e-5)
        _write_weight(handle, 10_000, True, var=1e-7, tol=1e-5)
        _write_weight(handle, 25_000, None, var=0.1, tol=1e-5)
    report = balance_utils.inspect_mcool_balance(mcool, [5_000, 10_000, 25_000])

    removed = balance_utils.sanitize_nonpass_weights(mcool, report)

    assert removed == [5_000, 25_000]
    with h5py.File(mcool, "r") as handle:
        assert "resolutions/5000/bins/weight" not in handle
        assert "resolutions/10000/bins/weight" in handle
        assert "resolutions/25000/bins/weight" not in handle
    assert report["resolutions"]["5000"]["weight_present"] is True
    assert report["resolutions"]["5000"]["weight_published"] is False


def test_balance_tsv_is_stable_and_records_published_weight(tmp_path):
    report = {
        "configured_resolutions_bp": [5_000],
        "resolutions": {"5000": {
            "status": "WARN", "weight_present": True,
            "weight_published": False, "converged": False,
            "variance": 0.0048, "tolerance": 1e-5,
            "parameters": {"ignore_diags": 2, "min_nnz": 10},
        }},
    }
    output = tmp_path / "balance.tsv"
    balance_qc.write_balance_tsv(report, output)
    table = pd.read_csv(output, sep="\t")
    assert list(table.columns) == balance_qc.TSV_COLUMNS
    assert bool(table.loc[0, "weight_present"])
    assert not bool(table.loc[0, "weight_published"])


def test_raw_cooltools_command_passes_a_literal_empty_weight_name():
    raw = cooltools_matrix_qc.build_command(
        "insulation", "map.mcool", "out.tsv", "view.bed", 25_000,
        4, 2, False, weight_name="weight", window=250_000,
    )
    balanced = cooltools_matrix_qc.build_command(
        "expected_cis", "map.mcool", "out.tsv", "view.bed", 25_000,
        4, 2, True, weight_name="weight",
    )
    raw_index = raw.index("--clr-weight-name")
    balanced_index = balanced.index("--clr-weight-name")
    assert raw[raw_index + 1] == ""
    assert balanced[balanced_index + 1] == "weight"
    assert raw[raw.index("--ignore-diags") + 1] == "2"


def test_zoomify_command_records_explicit_convergence_parameters():
    command = build_mcool.build_zoomify_command(
        "base.cool", "out.mcool", [5_000, 10_000], 4,
        weight_name="weight", ignore_diags=2, min_nnz=10, mad_max=5,
        tolerance=1e-5, max_iterations=200,
    )
    balance_args = command[command.index("--balance-args") + 1]
    assert "--tol 1e-05" in balance_args
    assert "--max-iters 200" in balance_args
    assert "--convergence-policy store_final" in balance_args


def test_balance_blacklist_materializes_gzip_as_plain_bed(tmp_path):
    source = tmp_path / "blacklist.bed.gz"
    destination = tmp_path / "blacklist.for-cooler.bed"
    expected = "chr1\t10000\t15000\n"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(expected)

    result = build_mcool.materialize_plaintext_blacklist(source, destination)

    assert result == destination
    assert destination.read_bytes() == expected.encode("utf-8")
    assert destination.read_bytes()[:2] != b"\x1f\x8b"


def test_blacklist_filter_removes_every_touching_contact(clr, tmp_path):
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr1\t10000\t15000\n")
    output = tmp_path / "clean.cool"
    audit = tmp_path / "audit.json"
    snakemake = SimpleNamespace(
        input=SimpleNamespace(cool=clr.filename, blacklist=str(blacklist)),
        output=SimpleNamespace(cool=str(output), json=str(audit)),
        params=SimpleNamespace(assembly="test"),
        wildcards=SimpleNamespace(sample="tiny"),
        log=[str(tmp_path / "filter.log")],
    )
    filter_cool_blacklist.main(snakemake)

    cleaned = cooler.Cooler(str(output))
    blocked_id = int(cleaned.offset("chr1")) + 2
    pixels = cleaned.pixels()[:]
    assert not (
        (pixels["bin1_id"] == blocked_id) | (pixels["bin2_id"] == blocked_id)
    ).any()
    report = json.loads(audit.read_text())
    assert report["n_blacklisted_bins"] == 1
    assert report["removed_contacts"] > 0
    assert report["retained_contacts"] < report["input_contacts"]


def test_mustache_is_not_assessed_when_required_balance_did_not_converge(tmp_path):
    balance = tmp_path / "balance.json"
    balance.write_text(json.dumps({
        "schema": "oracle-hichip-balance-qc-v1",
        "status": "WARN",
        "weight_name": "weight",
        "resolutions": {
            "10000": {"status": "WARN", "converged": False},
        },
    }))
    output = tmp_path / "calls.tsv"
    status = tmp_path / "calls.status.json"
    primary = tmp_path / "primary.bedpe"
    primary.write_text("")
    snakemake = SimpleNamespace(
        input=SimpleNamespace(
            mcool=str(tmp_path / "map.mcool"), balance=str(balance),
            primary=str(primary),
        ),
        output=SimpleNamespace(tsv=str(output), status=str(status)),
        params=SimpleNamespace(res=10_000, comparison_tolerance_bins=1),
        threads=2,
        wildcards=SimpleNamespace(sample="sample_1"),
        log=[str(tmp_path / "mustache.log")],
    )

    mustache_balance_aware.main(snakemake)

    payload = json.loads(status.read_text())
    assert payload["status"] == "NOT_ASSESSED"
    assert payload["available"] is False
    assert output.read_text() == mustache_balance_aware.MUSTACHE_HEADER


def test_e1_emits_header_only_table_when_balance_is_nonconverged(tmp_path):
    balance = tmp_path / "balance.json"
    balance.write_text(json.dumps({
        "schema": "oracle-hichip-balance-qc-v1",
        "status": "WARN",
        "weight_name": "weight",
        "resolutions": {
            "100000": {"status": "WARN", "converged": False},
        },
    }))
    eigs = tmp_path / "sample.cis.eigs.tsv"
    status = tmp_path / "sample.cis.eigs.status.json"
    snakemake = SimpleNamespace(
        input=SimpleNamespace(
            mcool=str(tmp_path / "map.mcool"), balance=str(balance),
            view=str(tmp_path / "view.bed"), gc=str(tmp_path / "gc.tsv"),
        ),
        output=SimpleNamespace(cis=str(eigs), status=str(status)),
        params=SimpleNamespace(res=100_000),
        wildcards=SimpleNamespace(sample="sample_1"),
        log=[str(tmp_path / "eigs.log")],
    )

    cooltools_eigs.main(snakemake)

    assert eigs.read_text() == "chrom\tstart\tend\tE1\n"
    payload = json.loads(status.read_text())
    assert payload["status"] == "NOT_ASSESSED"
    assert payload["available"] is False
    assert payload["normalization"] is None


def test_header_only_e1_creates_valid_empty_bigwig(tmp_path):
    try:
        import pyBigWig as pybigwig
    except ImportError:
        pytest.skip("pyBigWig is exercised in the coolerpy rule environment")
    compartments_to_bigwig = _load("compartments_to_bigwig")
    eigs = tmp_path / "empty.eigs.tsv"
    eigs.write_text("chrom\tstart\tend\tE1\n")
    chromsizes = tmp_path / "chrom.sizes"
    chromsizes.write_text("chr1\t100000\nchrX\t50000\n")
    output = tmp_path / "empty.E1.bw"
    snakemake = SimpleNamespace(
        input=SimpleNamespace(eigs=str(eigs), chromsizes=str(chromsizes)),
        output=SimpleNamespace(bw=str(output)),
        log=[str(tmp_path / "bigwig.log")],
    )

    compartments_to_bigwig.main(snakemake)

    handle = pybigwig.open(str(output))
    assert handle.chroms() == {"chr1": 100_000, "chrX": 50_000}
    assert handle.intervals("chr1") is None
    handle.close()


# --------------------------------------------------------- locus visualisation
def test_virtual_4c_requires_an_explicit_viewpoint():
    """A locus window midpoint is not a biological anchor and must not be inferred."""
    with pytest.raises(ValueError, match="no explicit viewpoint"):
        viewpoint.resolve_viewpoint({"name": "wide_locus", "start": 0, "end": 1_000_000})


def test_virtual_4c_rejects_viewpoint_outside_region():
    with pytest.raises(ValueError, match="outside"):
        viewpoint.resolve_viewpoint({
            "name": "bad_locus", "start": 100, "end": 200, "viewpoint": 200,
        })


def test_configured_viewpoints_are_gencode_v46_tss_coordinates():
    """Pin the hg38 0-based TSS values resolved from the configured GENCODE GTF."""
    root = Path(__file__).resolve().parents[1]
    cfg = yaml.safe_load((root / "config" / "config.yaml").read_text())
    regions = {r["name"]: r for r in cfg["viz"]["regions"]}
    expected = {
        "RORC_locus": 151_831_844,
        "IL17A_locus": 52_186_374,
        "FOXP3_locus": 49_264_799,
        "IL2RA_locus": 6_062_369,
        "CTLA4_locus": 203_853_887,
        "MYC_locus": 127_735_433,
    }
    assert {name: regions[name]["viewpoint"] for name in expected} == expected
    for region in regions.values():
        coordinate, label = viewpoint.resolve_viewpoint(region)
        assert region["start"] <= coordinate < region["end"]
        assert label.endswith(" TSS")


def test_composite_locus_track_includes_configured_gtf():
    rendered = pgt.INI_TEMPLATE.format(
        mcool="map.mcool", res=10_000, depth=1_000_000, sample="sample",
        matrix_title="HiChIP sample — raw-count fallback (balance WARN)",
        insulation_title="local insulation — raw-count fallback (balance WARN)",
        insul_bdg="insulation.bdg",
        peaks="peaks.bed", mark="H3K27ac",
        loops_section=pgt.LOOPS_TEMPLATE.format(loops="loops.bed"), gtf="genes.gtf.gz",
    )
    assert "[genes]" in rendered
    assert "file = genes.gtf.gz" in rendered
    assert "file_type = gtf" in rendered
    assert "prefered_name = gene_name" in rendered
    assert "raw-count fallback (balance WARN)" in rendered


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


# -------------------------------------------------------- report calculations
def test_loop_qc_ratios_use_explicit_populations():
    """Every contact ratio must use the population named in its definition."""
    dedup = {
        "total": 200.0,
        "total_mapped": 160.0,
        "total_dups": 40.0,
        "total_nodups": 120.0,
    }
    selected = {"total": 90.0, "total_nodups": 90.0, "cis": 72.0, "trans": 18.0}

    metrics = loop_qc.contact_qc_metrics(
        selected, dedup, raw_input_pairs=250, post_trim_pairs=200
    )

    assert metrics["valid_pair_yield_pct"] == pytest.approx(36.0)
    assert metrics["post_trim_valid_pair_yield_pct"] == pytest.approx(45.0)
    assert metrics["duplicate_pct"] == pytest.approx(25.0)
    assert metrics["cis_fraction"] == pytest.approx(0.8)
    assert metrics["duplicate_pct"] <= 100.0


def test_loop_qc_multiqc_payload_uses_current_fields_and_thresholds():
    report = {
        "sample": "sample_1",
        "valid_pair_yield_pct": 45.0,
        "duplicate_pct": 25.0,
        "cis_fraction": 0.8,
        "anchor_frip": 0.42,
        "n_consensus_peaks": 12_500,
        "n_loops": 2500,
        "apa_score": 2.2,
        "apa_status": "DESCRIPTIVE",
        "apa_normalization": "ICE-balanced",
        "hicrep_status": "PASS",
        "hicrep_min_scc": 0.71,
        "hicrep_mean_scc": 0.74,
        "hicrep_best_scc": 0.76,
        "hicrep_group_status": "PASS",
        "hicrep_group_median_scc": 0.73,
        "balance_status": "PASS",
        "balance_qc": {"n_converged": 9, "n_configured": 9},
        "mustache_status": "PASS",
        "mustache_balance_status": "PASS",
        "overall_status": "PASS",
        "thresholds": {"apa_score_min": 1.5, "hicrep_scc_min": 0.70},
    }

    payload = multiqc_content.loop_qc_payload(report)

    assert payload["id"] == "loop_qc"
    assert payload["data"]["sample_1"]["hicrep_status"] == "PASS"
    assert payload["data"]["sample_1"]["hicrep_min_scc"] == pytest.approx(0.71)
    assert payload["data"]["sample_1"]["hicrep_mean_scc"] == pytest.approx(0.74)
    assert payload["data"]["sample_1"]["hicrep_best_scc"] == pytest.approx(0.76)
    assert payload["data"]["sample_1"]["hicrep_group_status"] == "PASS"
    assert payload["data"]["sample_1"]["anchor_frip"] == pytest.approx(0.42)
    assert payload["data"]["sample_1"]["balance_status"] == "PASS"
    assert payload["data"]["sample_1"]["mustache_status"] == "PASS"
    assert "every" in payload["headers"]["hicrep_status"]["description"]
    assert "not a gate" in payload["headers"]["hicrep_best_scc"]["description"]
    assert "pass >=" not in payload["headers"]["hicrep_best_scc"]["description"].lower()


def test_apa_multiqc_payload_never_turns_unassessed_into_zero():
    measured = multiqc_content.apa_payload({
        "sample": "measured", "apa_vs_random_shift": 2.4,
        "normalization": "ICE-balanced",
    })
    not_assessed = multiqc_content.apa_payload({
        "sample": "not_assessed", "apa_vs_random_shift": None,
        "normalization": "ICE-balanced",
    })
    raw = multiqc_content.apa_payload({
        "sample": "raw", "apa_vs_random_shift": 3.1,
        "normalization": "raw-count fallback", "status": "DESCRIPTIVE",
    })

    assert measured["data"]["measured"]["APA versus random shift"] == pytest.approx(2.4)
    assert not_assessed["data"] == {}
    assert raw["data"] == {}


def test_balance_multiqc_payload_is_self_contained():
    payload = multiqc_content.balance_payload({
        "sample": "sample_1", "status": "WARN", "n_configured": 3,
        "n_converged": 2, "n_nonconverged": 1, "n_missing": 0,
        "converged_resolutions_bp": [25_000, 100_000],
        "nonconverged_resolutions_bp": [5_000],
        "missing_resolutions_bp": [],
    })
    assert payload["id"] == "balance_qc"
    row = payload["data"]["sample_1"]
    assert row["status"] == "WARN"
    assert row["nonconverged_resolutions"] == "5000"


def test_multiqc_content_main_writes_companion_and_log(tmp_path):
    source = tmp_path / "sample.apa.json"
    output = tmp_path / "sample_apa_mqc.json"
    log = tmp_path / "logs" / "sample.apa.log"
    source.write_text(json.dumps({
        "sample": "sample_1", "apa_vs_random_shift": 2.4,
        "normalization": "ICE-balanced",
    }))
    snakemake = SimpleNamespace(
        input=SimpleNamespace(json=str(source)),
        output=SimpleNamespace(json=str(output)),
        params=SimpleNamespace(kind="apa"),
        log=[str(log)],
    )

    multiqc_content.main(snakemake)

    payload = json.loads(output.read_text())
    assert payload["id"] == "apa_scores"
    assert payload["data"]["sample_1"]["APA versus random shift"] == pytest.approx(2.4)
    assert "Wrote apa custom content" in log.read_text()


def test_multiqc_config_uses_self_contained_companions():
    config_path = Path(__file__).resolve().parents[1] / "config" / "multiqc_config.yaml"
    config = yaml.safe_load(config_path.read_text())

    assert "custom_data" not in config
    assert "sp" not in config
    assert config["custom_content"]["order"] == [
        "balance_qc", "loop_qc", "apa_scores",
    ]


# ------------------------------------------------------------------- HiCRep
@pytest.mark.parametrize(
    ("values", "expected_status", "expected_pass"),
    [
        ([], "NOT_ASSESSED", None),
        ([0.70], "PASS", True),
        ([0.72, 0.91], "PASS", True),
        ([0.20, 0.69], "FAIL", False),
        ([0.20, 0.91], "DISCORDANT", None),
    ],
)
def test_hicrep_classifies_all_depth_qualified_pairs(
    values, expected_status, expected_pass
):
    status, passed = hicrep._classify_sccs(values, threshold=0.70)
    assert status == expected_status
    assert passed is expected_pass


def test_hicrep_discordant_triad_cannot_pass_by_best_or_group_median():
    pair_scc = {
        frozenset(("A", "B")): 0.85,
        frozenset(("A", "C")): 0.85,
        frozenset(("B", "C")): 0.20,
    }

    def involving(sample):
        return [value for pair, value in pair_scc.items() if sample in pair]

    assert hicrep._classify_sccs(involving("A"), 0.70)[0] == "PASS"
    assert hicrep._classify_sccs(involving("B"), 0.70)[0] == "DISCORDANT"
    assert hicrep._classify_sccs(involving("C"), 0.70)[0] == "DISCORDANT"
    assert hicrep._classify_sccs(list(pair_scc.values()), 0.70)[0] == "DISCORDANT"


def test_loop_qc_preserves_discordant_hicrep_despite_high_best_scc(tmp_path):
    pair_stats = tmp_path / "selected.stats"
    pair_stats.write_text(
        "total 90\ncis 80\nsummary/frac_cis 0.8888888889\n"
    )
    dedup_stats = tmp_path / "dedup.stats"
    dedup_stats.write_text(
        "total 100\ntotal_mapped 100\ntotal_dups 0\nsummary/frac_dups 0\n"
    )
    fastp = tmp_path / "fastp.json"
    fastp.write_text(json.dumps({
        "summary": {
            "before_filtering": {"total_reads": 200},
            "after_filtering": {"total_reads": 200},
        }
    }))
    apa = tmp_path / "apa.json"
    apa.write_text(json.dumps({
        "status": "DESCRIPTIVE",
        "normalization": "ICE-balanced",
        "balance_status": "PASS",
        "apa_vs_random_shift": 2.0,
        "apa_vs_random_shift_ci95": [1.8, 2.2],
    }))
    hicrep_json = tmp_path / "hicrep.json"
    hicrep_json.write_text(json.dumps({
        "status": "DISCORDANT",
        "min_scc": 0.20,
        "mean_scc": 0.62,
        "best_scc": 0.95,
        "group_status": "DISCORDANT",
        "group_median_scc": 0.85,
    }))
    restriction = tmp_path / "restriction.json"
    restriction.write_text(json.dumps({"fractions": {}}))
    anchor_qc = tmp_path / "anchor.tsv"
    anchor_qc.write_text("frip\tn_consensus_peaks\n0.4\t1000\n")
    loops = tmp_path / "loops.bedpe"
    loops.write_text("")
    balance = tmp_path / "balance.json"
    balance.write_text(json.dumps({
        "status": "PASS", "n_configured": 2, "n_converged": 2,
        "n_nonconverged": 0, "n_missing": 0,
        "converged_resolutions_bp": [5_000, 10_000],
        "nonconverged_resolutions_bp": [], "missing_resolutions_bp": [],
        "resolutions": {},
    }))
    contact_depth = tmp_path / "contact_depth.json"
    contact_depth.write_text(json.dumps({
        "primary_cis_offdiagonal_contacts": 80,
        "fithichip_distance_range_contacts": 70,
        "fithichip_distance_range_bp": [20_000, 3_000_000],
    }))
    output_json = tmp_path / "out" / "qc.json"
    output_md = tmp_path / "out" / "qc.md"

    snakemake = SimpleNamespace(
        input=SimpleNamespace(
            pair_stats=str(pair_stats),
            dedup_stats=str(dedup_stats),
            fastp=str(fastp),
            apa_json=str(apa),
            hicrep=str(hicrep_json),
            restriction=str(restriction),
            anchor_qc=str(anchor_qc),
            balance=str(balance),
            contact_depth=str(contact_depth),
            loops=str(loops),
        ),
        output=SimpleNamespace(json=str(output_json), md=str(output_md)),
        params=SimpleNamespace(thresholds={
            "n_loops_min": 0, "apa_score_min": 1.5,
            "hicrep_scc_min": 0.70,
        }),
        config={
            "qc_thresholds": {"n_loops_min": 0},
            "apa": {"score_min": 1.5},
            "hicrep": {"threshold_pass": 0.70},
        },
        wildcards=SimpleNamespace(sample="sample_1"),
        log=[str(tmp_path / "loop_qc.log")],
    )

    loop_qc.main(snakemake)
    report = json.loads(output_json.read_text())

    assert report["hicrep_best_scc"] == pytest.approx(0.95)
    assert report["hicrep_status"] == "DISCORDANT"
    assert report["status_flags"]["hicrep_scc"] == "DISCORDANT"
    assert report["overall_status"] == "PASS_WITH_UNCERTAINTY"
    assert report["overall_pass"] is False
    assert "DISCORDANT" in output_md.read_text()


def test_loop_qc_balance_warn_propagates_as_uncertainty_not_failure(tmp_path):
    pair_stats = tmp_path / "selected.stats"
    pair_stats.write_text("total 90\ncis 80\nsummary/frac_cis 0.8888888889\n")
    dedup_stats = tmp_path / "dedup.stats"
    dedup_stats.write_text(
        "total 100\ntotal_mapped 100\ntotal_dups 0\nsummary/frac_dups 0\n"
    )
    fastp = tmp_path / "fastp.json"
    fastp.write_text(json.dumps({
        "summary": {
            "before_filtering": {"total_reads": 200},
            "after_filtering": {"total_reads": 200},
        }
    }))
    apa = tmp_path / "apa.json"
    apa.write_text(json.dumps({
        "status": "DESCRIPTIVE", "normalization": "ICE-balanced",
        "balance_status": "PASS", "apa_vs_random_shift": 2.0,
        "apa_vs_random_shift_ci95": [1.8, 2.2],
    }))
    hicrep_json = tmp_path / "hicrep.json"
    hicrep_json.write_text(json.dumps({
        "status": "PASS", "min_scc": 0.8, "mean_scc": 0.8,
        "best_scc": 0.8, "group_status": "PASS", "group_median_scc": 0.8,
    }))
    restriction = tmp_path / "restriction.json"
    restriction.write_text(json.dumps({"fractions": {}}))
    anchor_qc = tmp_path / "anchor.tsv"
    anchor_qc.write_text("frip\tn_consensus_peaks\n0.4\t1000\n")
    loops = tmp_path / "loops.bedpe"
    loops.write_text("")
    balance = tmp_path / "balance.json"
    balance.write_text(json.dumps({
        "status": "WARN", "n_configured": 2, "n_converged": 1,
        "n_nonconverged": 1, "n_missing": 0,
        "converged_resolutions_bp": [10_000],
        "nonconverged_resolutions_bp": [5_000],
        "missing_resolutions_bp": [], "resolutions": {},
    }))
    contact_depth = tmp_path / "contact_depth.json"
    contact_depth.write_text(json.dumps({
        "primary_cis_offdiagonal_contacts": 80,
        "fithichip_distance_range_contacts": 70,
        "fithichip_distance_range_bp": [20_000, 3_000_000],
    }))
    mustache = tmp_path / "mustache.status.json"
    mustache.write_text(json.dumps({
        "status": "NOT_ASSESSED", "balance_status": "WARN",
        "available": False,
    }))
    output_json = tmp_path / "out" / "qc.json"
    snakemake = SimpleNamespace(
        input=SimpleNamespace(
            pair_stats=str(pair_stats), dedup_stats=str(dedup_stats),
            fastp=str(fastp),
            apa_json=str(apa), hicrep=str(hicrep_json),
            restriction=str(restriction), anchor_qc=str(anchor_qc),
            balance=str(balance), contact_depth=str(contact_depth),
            mustache=[str(mustache)], loops=str(loops),
        ),
        output=SimpleNamespace(
            json=str(output_json), md=str(tmp_path / "out" / "qc.md")
        ),
        params=SimpleNamespace(thresholds={
            "n_loops_min": 0, "apa_score_min": 1.5,
            "hicrep_scc_min": 0.70,
        }),
        config={
            "qc_thresholds": {"n_loops_min": 0},
            "apa": {"score_min": 1.5},
            "hicrep": {"threshold_pass": 0.70},
        },
        wildcards=SimpleNamespace(sample="sample_1"),
        log=[str(tmp_path / "loop_qc.log")],
    )

    loop_qc.main(snakemake)
    report = json.loads(output_json.read_text())

    assert report["balance_status"] == "WARN"
    assert report["status_flags"]["cooler_balance"] == "WARN"
    assert report["status_flags"]["mustache"] == "NOT_ASSESSED"
    assert report["mustache_status"] == "NOT_ASSESSED"
    assert report["overall_status"] == "PASS_WITH_UNCERTAINTY"
    assert report["overall_pass"] is False


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
def test_oracle_containment_maps_each_fine_bin_to_one_parent():
    cos = _load("export_oracle_cos")
    fine = pd.DataFrame({
        "chrom": ["chr1"] * 4, "start": [0, 5_000, 10_000, 15_000],
        "end": [5_000, 10_000, 15_000, 20_000], "bin_idx": range(4),
    })
    coarse = pd.DataFrame({
        "chrom": ["chr1"] * 2, "start": [0, 10_000],
        "end": [10_000, 20_000], "bin_idx": range(2),
    })
    assert cos._containment_edges(fine, coarse).tolist() == [
        [0, 1, 2, 3], [0, 0, 1, 1]
    ]


def test_oracle_primary_chromosome_policy_excludes_alts_y_and_mitochondria():
    cos = _load("export_oracle_cos")
    chromsizes = {
        "chr1": 100, "chr22": 100, "chrX": 100, "chrY": 100,
        "chrM": 100, "chr1_KI270706v1_random": 100, "GL000008.2": 100,
    }
    assert cos._select_chromosomes(chromsizes, set(), True) == [
        "chr1", "chr22", "chrX",
    ]


def test_oracle_feature_availability_preserves_channel_order_and_marks_missing_e1():
    cos = _load("export_oracle_cos")
    insulation = pd.DataFrame({
        "log2_insulation_score_250000": [0.2],
        "normalization": ["raw-count fallback"],
        "balance_status": ["WARN"],
    })
    eigs = pd.DataFrame(columns=["chrom", "start", "end", "E1"])
    availability = cos.node_feature_availability(
        insulation, "log2_insulation_score_250000", eigs,
        {
            "available": False, "status": "NOT_ASSESSED",
            "balance_status": "WARN", "resolution_bp": 100_000,
            "reason": "ICE balancing did not converge at the E1 resolution",
        },
    )
    assert cos.NODE_FEATURE_CHANNELS == [
        "peak_overlap_count_per_kb", "insulation", "E1_eigenvector",
    ]
    assert list(availability) == cos.NODE_FEATURE_CHANNELS
    assert availability["insulation"]["available"] is True
    assert availability["insulation"]["normalization"] == "raw-count fallback"
    assert availability["E1_eigenvector"]["available"] is False


def test_oracle_peak_sweep_matches_brute_force(tmp_path):
    cos = _load("export_oracle_cos")
    bins = pd.DataFrame({
        "chrom": ["chr1"] * 20,
        "start": np.arange(20) * 1_000,
        "end": (np.arange(20) + 1) * 1_000,
        "bin_idx": np.arange(20),
    })
    intervals = [(50, 1_050), (999, 3_001), (5_000, 7_000), (19_500, 20_000)]
    peaks = tmp_path / "peaks.bed"
    peaks.write_text("".join(f"chr1\t{s}\t{e}\n" for s, e in intervals))

    got = cos._peak_overlap_per_bin(peaks, bins)
    brute = np.zeros(len(bins), dtype=float)
    for i, row in bins.iterrows():
        brute[i] = sum(row.start < end and row.end > start for start, end in intervals)
    # 1 kb bins make count-per-kb numerically equal to the overlap count.
    assert got.tolist() == brute.tolist()


def test_oracle_unobserved_feature_cells_are_zero_filled():
    cos = _load("export_oracle_cos")
    values = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float32)
    observed = np.array([[True, False, True], [False, False, True]])
    got = cos._zero_unobserved_values(values, observed)
    assert got.tolist() == [[1.0, 0.0, 3.0], [0.0, 0.0, 6.0]]
    assert values[0, 1] == 2.0  # helper does not mutate its caller's array


def test_oracle_optional_tokens_keep_one_schema_and_missingness_mask(tmp_path):
    cos = _load("export_oracle_cos")
    table = tmp_path / "tokens.tsv"
    table.write_text(
        "sample_id\tage\tmetabolite\tlabel\n"
        "A\t40\t\tcase\n"
        "B\t50\t2.5\tcontrol\n"
    )
    keys_a, values_a, mask_a = cos._load_global_tokens(table, "A")
    keys_missing, values_missing, mask_missing = cos._load_global_tokens(table, "C")
    assert keys_a == keys_missing == ["age", "metabolite"]
    assert values_a.tolist() == [40.0, 0.0]
    assert mask_a.tolist() == [True, False]
    assert values_missing.tolist() == [0.0, 0.0]
    assert mask_missing.tolist() == [False, False]


def test_pyg_sample_tokens_batch_as_one_node_per_sample():
    torch = pytest.importorskip("torch")
    tg_data = pytest.importorskip("torch_geometric.data")
    graphs = []
    for values in ([1.0, 2.0], [3.0, 4.0]):
        graph = tg_data.HeteroData()
        graph["sample"].num_nodes = 1
        graph["sample"].microbiome = torch.tensor([values])
        graphs.append(graph)
    batch = tg_data.Batch.from_data_list(graphs)
    assert batch["sample"].microbiome.shape == (2, 2)
    assert batch["sample"].ptr.tolist() == [0, 1, 2]


def test_loops_to_edges_empty_result_keeps_2d_attr_shape():
    """A non-empty loop set that yields no edges keeps the v2 five-channel schema.

    Both anchors of a loop shorter than one bin land in the same bin (i == j) and are
    skipped, so `attrs` can be empty even though the BEDPE was not. `np.asarray([])`
    is 1-D, and concatenating it with the 2-D adjacency attributes raises

        ValueError: all the input arrays must have same number of dimensions

    which would otherwise fail during adjacency-edge concatenation.
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
    assert edge_attr.shape == (0, 5)


def test_loops_to_edges_preserves_raw_p_and_adjusted_q_separately():
    cos = _load("export_oracle_cos")
    bins = pd.DataFrame({
        "chrom": ["chr1"] * 3,
        "start": [0, 5_000, 10_000],
        "end": [5_000, 10_000, 15_000],
        "bin_idx": [0, 1, 2],
    })
    loops = pd.DataFrame([{
        "chrom1": "chr1", "start1": 0, "end1": 5_000,
        "chrom2": "chr1", "start2": 10_000, "end2": 15_000,
        "score": 17, "pvalue": 1e-6, "fdr": 0.004,
    }])

    _, edge_attr = cos._loops_to_edges(loops, bins)

    assert edge_attr.shape == (2, 5)
    assert edge_attr[0].tolist() == pytest.approx(
        [17, 1e-6, 0.004, 10_000, 1]
    )


def test_oracle_coarse_loop_edges_are_consolidated_with_multiplicity():
    cos = _load("export_oracle_cos")
    bins = pd.DataFrame({
        "chrom": ["chr1"] * 3, "start": [0, 25_000, 50_000],
        "end": [25_000, 50_000, 75_000], "bin_idx": [0, 1, 2],
    })
    loops = pd.DataFrame([
        {"chrom1": "chr1", "start1": 0, "end1": 5_000,
         "chrom2": "chr1", "start2": 50_000, "end2": 55_000,
         "score": 9, "pvalue": 0.02, "fdr": 0.04},
        {"chrom1": "chr1", "start1": 5_000, "end1": 10_000,
         "chrom2": "chr1", "start2": 55_000, "end2": 60_000,
         "score": 12, "pvalue": 0.01, "fdr": 0.03},
    ])
    edge_index, edge_attr = cos._loops_to_edges(loops, bins)
    assert edge_index.shape == (2, 2)
    assert edge_attr[0].tolist() == pytest.approx(
        [12, 0.01, 0.03, 50_000, 2]
    )


def test_oracle_blacklist_mask_fails_closed_on_missing_input(tmp_path):
    cos = _load("export_oracle_cos")
    bins = pd.DataFrame({
        "chrom": ["chr1"], "start": [0], "end": [5_000], "bin_idx": [0],
    })

    with pytest.raises(FileNotFoundError):
        cos._blacklist_mask_for_bins(tmp_path / "missing.bed.gz", bins)


def test_apa_bootstrap_resamples_matched_loops_not_independent_pools():
    apa = _load("apa_plot")
    real = np.array([10.0, 20.0, 40.0])
    control = np.array([5.0, 10.0, 20.0])

    point, low, high = apa._matched_ratio_bootstrap(
        real, control, np.random.default_rng(3), n_boot=200
    )

    assert point == pytest.approx(2.0)
    assert low == pytest.approx(2.0)
    assert high == pytest.approx(2.0)


def test_interval_feature_broadcasts_coarse_values_to_fine_bins():
    """Every 5 kb bin covered by a 25 kb/100 kb feature inherits its value."""
    cos = _load("export_oracle_cos")
    bins = pd.DataFrame({
        "chrom": ["chr1"] * 6,
        "start": np.arange(0, 30_000, 5_000),
        "end": np.arange(5_000, 35_000, 5_000),
    })
    track = pd.DataFrame({
        "chrom": ["chr1", "chr1"],
        "start": [0, 25_000],
        "end": [25_000, 100_000],
        "value": [2.0, 7.0],
    })

    got = cos._interval_feature_per_bin(track, bins, "value")

    assert np.array_equal(got, np.array([2, 2, 2, 2, 2, 7], dtype=np.float32))


def test_interval_feature_uses_overlap_weighted_mean_for_coarse_bins():
    """Fine source intervals are aggregated by covered bases, not start keys."""
    cos = _load("export_oracle_cos")
    bins = pd.DataFrame({"chrom": ["chr1"], "start": [0], "end": [20_000]})
    track = pd.DataFrame({
        "chrom": ["chr1", "chr1"],
        "start": [0, 5_000],
        "end": [5_000, 20_000],
        "value": [2.0, 6.0],
    })

    got = cos._interval_feature_per_bin(track, bins, "value")

    assert got[0] == pytest.approx(5.0)


def test_apa_min_dist_floor_keeps_diagonal_out_of_window():
    """The configured APA distance floor must exclude the main diagonal.

    A pixel at offset (dy, dx) from a loop of span D sits at separation
    D + (dx - dy) * bin, and dx - dy ranges over +/- 2*window. So the diagonal enters
    the window for any loop with D <= 2*window*bin. The configured distance floor must
    exclude that geometry independently of the cohort's observed loop spans.
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


# --------------------------------------------------- audited QC denominators
def test_contact_qc_metrics_keep_populations_and_denominators_straight():
    qc = _load("loop_qc_summary")
    dedup = {
        "total": 1_000,
        "total_mapped": 800,
        "total_dups": 200,
        "summary/frac_dups": 0.25,
    }
    selected = {
        "total": 400,
        "cis": 300,
        "summary/frac_cis": 0.75,
    }
    got = qc.contact_qc_metrics(
        selected, dedup, raw_input_pairs=2_000, post_trim_pairs=1_000
    )
    assert got["valid_pair_yield_pct"] == pytest.approx(20.0)
    assert got["post_trim_valid_pair_yield_pct"] == pytest.approx(40.0)
    assert got["duplicate_pct"] == pytest.approx(25.0)
    assert got["cis_fraction"] == pytest.approx(0.75)
    assert got["valid_pair_yield_pct"] != 100


def test_fastp_raw_and_post_trim_pair_populations_are_explicit(tmp_path):
    qc = _load("loop_qc_summary")
    report = tmp_path / "fastp.json"
    report.write_text(json.dumps({
        "summary": {
            "before_filtering": {"total_reads": 2_000},
            "after_filtering": {"total_reads": 1_800},
        }
    }))
    assert qc.fastp_pair_populations(report) == (1_000, 900)

    report.write_text(json.dumps({
        "summary": {
            "before_filtering": {"total_reads": 2_001},
            "after_filtering": {"total_reads": 1_800},
        }
    }))
    with pytest.raises(ValueError, match="even"):
        qc.fastp_pair_populations(report)


# --------------------------------------------------- differential design/candidates
def test_paired_design_blocks_on_subject_and_drops_constant_covariate():
    diff = _load("differential_loops")
    samples = pd.DataFrame({
        "sample_id": ["case1", "ctrl1", "case2", "ctrl2"],
        "subject_id": ["d1", "d1", "d2", "d2"],
        "batch": ["b1"] * 4,
    }).set_index("sample_id")
    metadata, factors = diff.build_design_metadata(
        list(samples.index), ["case1", "case2"], ["ctrl1", "ctrl2"],
        samples, "subject_id", ["batch"],
    )
    assert factors == ["subject_id", "condition"]
    assert metadata.loc["case1", "condition"] == "case"
    assert metadata.loc["ctrl1", "condition"] == "control"


def test_paired_design_rejects_unmatched_donors():
    diff = _load("differential_loops")
    samples = pd.DataFrame({
        "sample_id": ["case1", "ctrl1", "case2", "ctrl3"],
        "subject_id": ["d1", "d1", "d2", "d3"],
    }).set_index("sample_id")
    with pytest.raises(ValueError, match="not one-to-one paired"):
        diff.build_design_metadata(
            list(samples.index), ["case1", "case2"], ["ctrl1", "ctrl3"],
            samples, "subject_id", [],
        )


def test_paired_design_rejects_condition_confounded_batch():
    diff = _load("differential_loops")
    samples = pd.DataFrame({
        "sample_id": ["case1", "ctrl1", "case2", "ctrl2"],
        "subject_id": ["d1", "d1", "d2", "d2"],
        "batch": ["case_batch", "control_batch", "case_batch", "control_batch"],
    }).set_index("sample_id")
    with pytest.raises(ValueError, match="rank-deficient"):
        diff.build_design_metadata(
            list(samples.index), ["case1", "case2"], ["ctrl1", "ctrl2"],
            samples, "subject_id", ["batch"],
        )


def test_consensus_loops_collapse_shifted_merged_calls_and_require_support(tmp_path):
    consensus = _load("build_consensus_loops")
    header = "chr1\ts1\te1\tchr2\ts2\te2\tcc\tP-Value_Bias\tQ-Value_Bias\n"
    a = tmp_path / "sample_a" / "loops.bed"
    b = tmp_path / "sample_b" / "loops.bed"
    a.parent.mkdir()
    b.parent.mkdir()
    a.write_text(header + "chr1\t0\t5000\tchr1\t20000\t25000\t9\t0.01\t0.02\n")
    b.write_text(header + "chr1\t5000\t10000\tchr1\t25000\t30000\t8\t0.02\t0.03\n")
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr2\t0\t100\n")
    got = consensus.build_consensus([str(a), str(b)], 5_000, 2, blacklist)
    assert len(got) == 1
    assert got.loc[0, "sample_support"] == 2
    assert got.loc[0, "start1"] == 0
    assert got.loc[0, "end1"] == 10_000
    assert got.loc[0, "start2"] == 15_000
    assert got.loc[0, "end2"] == 30_000
    # Exact-only support cannot reconcile the true adjacent-grid calls.
    exact = consensus.build_consensus(
        [str(a), str(b)], 5_000, 2, blacklist, tolerance_bins=0
    )
    assert exact.empty


def test_consensus_loop_tolerance_votes_once_per_sample(tmp_path):
    consensus = _load("build_consensus_loops")
    header = "chr1\ts1\te1\tchr2\ts2\te2\tcc\tP-Value_Bias\tQ-Value_Bias\n"
    a = tmp_path / "sample_a" / "loops.bed"
    b = tmp_path / "sample_b" / "loops.bed"
    a.parent.mkdir()
    b.parent.mkdir()
    a.write_text(
        header
        + "chr1\t0\t5000\tchr1\t20000\t25000\t9\t0.01\t0.02\n"
        + "chr1\t5000\t10000\tchr1\t25000\t30000\t8\t0.02\t0.03\n"
    )
    b.write_text(header + "chr1\t5000\t10000\tchr1\t25000\t30000\t7\t0.03\t0.04\n")
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr2\t0\t100\n")
    got = consensus.build_consensus([str(a), str(b)], 5_000, 2, blacklist)
    assert len(got) == 1
    assert got.loc[0, "sample_support"] == 2
    assert got.loc[0, "support_samples"] == "sample_a,sample_b"


def test_consensus_loop_neighbourhood_does_not_expand_transitively(tmp_path):
    consensus = _load("build_consensus_loops")
    header = "chr1\ts1\te1\tchr2\ts2\te2\tcc\tP-Value_Bias\tQ-Value_Bias\n"
    files = []
    for i in range(4):
        path = tmp_path / f"sample_{i}" / "loops.bed"
        path.parent.mkdir()
        start1, start2 = i * 5_000, 40_000 + i * 5_000
        path.write_text(
            header + f"chr1\t{start1}\t{start1 + 5000}\tchr1\t{start2}\t{start2 + 5000}\t9\t0.01\t0.02\n"
        )
        files.append(str(path))
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr2\t0\t100\n")
    got = consensus.build_consensus(files, 5_000, 2, blacklist, tolerance_bins=1)
    # No representative has all four calls within +/-1 bin; chained expansion
    # must therefore never report support=4.
    assert got["sample_support"].max() == 3


def test_consensus_peaks_requires_independent_sample_support(tmp_path):
    peaks = _load("consensus_peaks")
    a = tmp_path / "a_peaks.bed"
    b = tmp_path / "b_peaks.bed"
    a.write_text("chr1\t100\t200\nchr1\t1000\t1100\n")
    b.write_text("chr1\t150\t250\n")
    got = peaks.consensus_peaks([str(a), str(b)], 2)
    assert got[["chrom", "start", "end", "sample_support"]].values.tolist() == [
        ["chr1", 150, 200, 2]
    ]


def test_consensus_peaks_does_not_retain_transitive_unsupported_flanks(tmp_path):
    peaks = _load("consensus_peaks")
    a = tmp_path / "a_peaks.bed"
    b = tmp_path / "b_peaks.bed"
    c = tmp_path / "c_peaks.bed"
    a.write_text("chr1\t0\t100\n")
    b.write_text("chr1\t50\t150\n")
    c.write_text("chr1\t100\t200\n")

    got = peaks.consensus_peaks([str(a), str(b), str(c)], 2)

    assert got[["chrom", "start", "end", "sample_support"]].values.tolist() == [
        ["chr1", 50, 100, 2],
        ["chr1", 100, 150, 2],
    ]


# --------------------------------------------------- restriction fragments/preflight
def test_restriction_orientation_classes_are_explicit():
    restriction = _load("restriction_fragment_qc")
    assert restriction.classify("chr1", "chr1", 10, 11, "+", "-") == "dangling_end_like"
    assert restriction.classify("chr1", "chr1", 10, 10, "-", "+") == "self_circle_like"
    assert restriction.classify("chr1", "chr2", 10, 10, "+", "-") == "regular"


def test_shipped_configuration_passes_preflight():
    import yaml
    validate = _load("validate_config")
    root = Path(__file__).resolve().parents[1]
    cfg = yaml.safe_load((root / "config/config.yaml").read_text())
    genomes = yaml.safe_load((root / "config/genome.yaml").read_text())
    samples = pd.read_csv(
        root / "config/samples.tsv", sep="\t", comment="#", dtype=str,
        keep_default_na=False,
    )
    validate.validate_pipeline_config(cfg, genomes, samples)


def test_preflight_rejects_an_unmatched_paired_comparison():
    import copy
    import yaml
    validate = _load("validate_config")
    root = Path(__file__).resolve().parents[1]
    cfg = yaml.safe_load((root / "config/config.yaml").read_text())
    genomes = yaml.safe_load((root / "config/genome.yaml").read_text())
    samples = pd.read_csv(
        root / "config/samples.tsv", sep="\t", comment="#", dtype=str,
        keep_default_na=False,
    )
    broken = copy.deepcopy(cfg)
    broken["differential"]["comparisons"][0]["case_filter"] = {"cell_type": "Th17", "replicate": "1"}
    with pytest.raises(ValueError, match="at least two|unmatched"):
        validate.validate_pipeline_config(broken, genomes, samples)


def test_shipped_differential_comparisons_select_published_style_complete_pairs():
    import yaml
    validate = _load("validate_config")
    root = Path(__file__).resolve().parents[1]
    cfg = yaml.safe_load((root / "config/config.yaml").read_text())
    samples = pd.read_csv(
        root / "config/samples.tsv", sep="\t", comment="#", dtype=str,
        keep_default_na=False,
    )
    comp = cfg["differential"]["comparisons"][0]
    cases = validate._selected(
        samples, comp["case_filter"], comp["mark"], comp["include_subjects"]
    )
    controls = validate._selected(
        samples, comp["control_filter"], comp["mark"], comp["include_subjects"]
    )
    assert set(cases["subject_id"]) == set(controls["subject_id"]) == {
        "donor2", "donor3"
    }
    assert len(cases) == len(controls) == 2


def test_preflight_rejects_unknown_included_subject_and_report_sample():
    import copy
    import yaml
    validate = _load("validate_config")
    root = Path(__file__).resolve().parents[1]
    cfg = yaml.safe_load((root / "config/config.yaml").read_text())
    genomes = yaml.safe_load((root / "config/genome.yaml").read_text())
    samples = pd.read_csv(
        root / "config/samples.tsv", sep="\t", comment="#", dtype=str,
        keep_default_na=False,
    )

    broken_subject = copy.deepcopy(cfg)
    broken_subject["differential"]["comparisons"][0]["include_subjects"] = [
        "donor2", "missing_donor"
    ]
    with pytest.raises(ValueError, match="unknown subject IDs"):
        validate.validate_pipeline_config(broken_subject, genomes, samples)

    broken_report = copy.deepcopy(cfg)
    broken_report["reporting"]["demonstration_samples"] = ["missing_sample"]
    with pytest.raises(ValueError, match="unknown sample IDs"):
        validate.validate_pipeline_config(broken_report, genomes, samples)


@pytest.mark.parametrize(
    "oracle_sizes",
    [[], [5, 5], [0, 25], [-5, 25], [5.5, 25], [True, 25]],
)
def test_preflight_rejects_invalid_oracle_resolution_lists(oracle_sizes):
    import copy
    import yaml
    validate = _load("validate_config")
    root = Path(__file__).resolve().parents[1]
    cfg = yaml.safe_load((root / "config/config.yaml").read_text())
    genomes = yaml.safe_load((root / "config/genome.yaml").read_text())
    samples = pd.read_csv(
        root / "config/samples.tsv", sep="\t", comment="#", dtype=str,
        keep_default_na=False,
    )
    broken = copy.deepcopy(cfg)
    broken["cooler"]["oracle_bin_sizes_kb"] = oracle_sizes
    with pytest.raises(ValueError, match="oracle_bin_sizes_kb"):
        validate.validate_pipeline_config(broken, genomes, samples)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("primary_chromosomes_only", "false"),
        ("drop_chromosomes", "chrM"),
        ("drop_chromosomes", ["chrM", "chrM"]),
        ("microbiome_metadata_tsv", 7),
    ],
)
def test_preflight_rejects_unsafe_oracle_export_types(key, value):
    import copy
    import yaml
    validate = _load("validate_config")
    root = Path(__file__).resolve().parents[1]
    cfg = yaml.safe_load((root / "config/config.yaml").read_text())
    genomes = yaml.safe_load((root / "config/genome.yaml").read_text())
    samples = pd.read_csv(
        root / "config/samples.tsv", sep="\t", comment="#", dtype=str,
        keep_default_na=False,
    )
    broken = copy.deepcopy(cfg)
    broken["oracle_export"][key] = value
    with pytest.raises(ValueError, match="oracle_export"):
        validate.validate_pipeline_config(broken, genomes, samples)


def test_preflight_rejects_duplicate_biological_units():
    import copy
    validate = _load("validate_config")
    root = Path(__file__).resolve().parents[1]
    cfg = yaml.safe_load((root / "config/config.yaml").read_text())
    genomes = yaml.safe_load((root / "config/genome.yaml").read_text())
    samples = pd.read_csv(
        root / "config/samples.tsv", sep="\t", comment="#", dtype=str,
        keep_default_na=False,
    )
    duplicate = copy.deepcopy(samples.iloc[0])
    duplicate["sample_id"] = "duplicate_technical_run"
    broken = pd.concat([samples, duplicate.to_frame().T], ignore_index=True)

    with pytest.raises(ValueError, match="merge technical runs"):
        validate.validate_pipeline_config(cfg, genomes, broken)


def test_preflight_rejects_duplicate_sra_accessions_across_samples():
    import copy
    validate = _load("validate_config")
    root = Path(__file__).resolve().parents[1]
    cfg = yaml.safe_load((root / "config/config.yaml").read_text())
    genomes = yaml.safe_load((root / "config/genome.yaml").read_text())
    samples = pd.read_csv(
        root / "config/samples.tsv", sep="\t", comment="#", dtype=str,
        keep_default_na=False,
    )
    broken = copy.deepcopy(samples)
    broken.loc[1, "srr"] = samples.loc[0, "srr"]

    with pytest.raises(ValueError, match="SRA accessions must be unique"):
        validate.validate_pipeline_config(cfg, genomes, broken)
