"""Regression tests for the q-free exact-grid differential contract."""
import copy
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))


def _load(name: str):
    spec = importlib.util.spec_from_file_location(
        f"differential_contract_{name}", SCRIPTS / f"{name}.py"
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


normalizer = _load("normalize_fithichip_all")
universe = _load("build_differential_universe")
diff = _load("differential_loops")
validate = _load("validate_config")
reporting_filter = _load("filter_fithichip_loops")


RAW_HEADER = (
    "chr1\ts1\te1\tchr2\ts2\te2\tcc\tisPeak1\tisPeak2\t"
    "P-Value_Bias\tQ-Value_Bias\n"
)
COORDS = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]


def _raw_row(
    chrom: str, start1: int, start2: int, count: int, qvalue: float
) -> str:
    return (
        f"{chrom}\t{start1}\t{start1 + 5000}\t{chrom}\t{start2}\t"
        f"{start2 + 5000}\t{count}\t1\t0\t0.2\t{qvalue}\n"
    )


def test_normalizer_preserves_high_q_rows_and_separates_abundance_filter(tmp_path):
    source = tmp_path / "sample.interactions_FitHiC.bed"
    source.write_text(
        RAW_HEADER
        + _raw_row("chr1", 0, 20_000, 4, 0.001)
        + _raw_row("chr1", 5_000, 30_000, 8, 0.90)
        + _raw_row("chr1", 10_000, 40_000, 9, 0.002)
        + (
            "chr1\t15000\t20000\tchr1\t50000\t55000\t7\t1\t0\t"
            "0.2\tNA\n"
        )
        + _raw_row("chrUn_GL0001", 0, 20_000, 7, 0.80)
    )
    chromsizes = tmp_path / "chrom.sizes"
    chromsizes.write_text("chr1\t100000\n")
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr1\t40000\t45000\n")
    all_output = tmp_path / "sample.all.tsv.gz"
    eligible_output = tmp_path / "sample.eligible.tsv.gz"
    audit_output = tmp_path / "sample.audit.json"

    audit = normalizer.normalize_all_interactions(
        source,
        blacklist,
        chromsizes,
        all_output,
        eligible_output,
        audit_output,
        sample="sample",
        bin_size=5_000,
        lower_distance=10_000,
        upper_distance=50_000,
        min_count=5,
        interaction_type="Peak-to-ALL",
        source_relative="native/sample.interactions_FitHiC.bed",
        chunk_size=2,
    )

    all_rows = pd.read_csv(all_output, sep="\t")
    eligible = pd.read_csv(eligible_output, sep="\t")
    assert all_rows["score"].tolist() == [4, 8, 7]
    assert eligible["score"].tolist() == [8, 7]
    assert eligible.loc[0, "fdr"] == pytest.approx(0.90)
    assert pd.isna(eligible.loc[1, "fdr"])
    assert audit["fithichip_q_filter"] is None
    assert audit["merge_nearby"] is False
    assert audit["retained_q_gt_0_05_rows"] == 1
    assert audit["retained_missing_q_rows"] == 1
    assert audit["removed_blacklist"] == 1
    assert audit["removed_non_primary"] == 1


def test_normalizer_rejects_non_native_grid(tmp_path):
    source = tmp_path / "bad.bed"
    source.write_text(RAW_HEADER + _raw_row("chr1", 1, 20_000, 8, 0.5))
    chromsizes = tmp_path / "chrom.sizes"
    chromsizes.write_text("chr1\t100000\n")
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("")
    with pytest.raises(ValueError, match="native bin"):
        normalizer.normalize_all_interactions(
            source,
            blacklist,
            chromsizes,
            tmp_path / "all.gz",
            tmp_path / "eligible.gz",
            tmp_path / "audit.json",
            sample="bad",
            bin_size=5_000,
            lower_distance=10_000,
            upper_distance=50_000,
            min_count=5,
            interaction_type="Peak-to-ALL",
            source_relative="bad.bed",
        )


@pytest.mark.parametrize("bad_value", ["garbage", "1.2", "inf"])
def test_normalizer_rejects_malformed_probability_values(tmp_path, bad_value):
    source = tmp_path / "bad_probability.bed"
    source.write_text(
        RAW_HEADER
        + f"chr1\t0\t5000\tchr1\t20000\t25000\t8\t1\t0\t0.2\t{bad_value}\n"
    )
    chromsizes = tmp_path / "chrom.sizes"
    chromsizes.write_text("chr1\t100000\n")
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("")
    with pytest.raises(ValueError, match="probabilities"):
        normalizer.normalize_all_interactions(
            source,
            blacklist,
            chromsizes,
            tmp_path / "all.gz",
            tmp_path / "eligible.gz",
            tmp_path / "audit.json",
            sample="bad_probability",
            bin_size=5_000,
            lower_distance=10_000,
            upper_distance=50_000,
            min_count=5,
            interaction_type="Peak-to-ALL",
            source_relative="bad_probability.bed",
        )


def _write_normalized(path: Path, rows: list[tuple]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows, columns=COORDS + ["score", "pvalue", "fdr"])
    frame.to_csv(path, sep="\t", index=False, compression="gzip")


def _universe_fixture(tmp_path: Path):
    samples = ["s1", "s2", "s3", "s4"]
    loop1 = ("chr1", 0, 5_000, "chr1", 20_000, 25_000)
    loop2 = ("chr1", 5_000, 10_000, "chr1", 30_000, 35_000)
    by_sample = {
        "s1": [(*loop1, 8, 0.2, 0.90), (*loop2, 1, 0.4, 0.95)],
        "s2": [(*loop1, 7, 0.3, 0.80)],
        "s3": [(*loop1, 3, 0.6, 0.99), (*loop2, 9, 0.2, 0.70)],
        "s4": [(*loop2, 10, 0.3, 0.60)],
    }
    all_files, eligible_files, audit_files = [], [], []
    for sample in samples:
        sample_dir = tmp_path / sample
        all_path = sample_dir / f"{sample}.interactions_FitHiC.all.tsv.gz"
        eligible_path = sample_dir / f"{sample}.interactions_FitHiC.eligible.tsv.gz"
        rows = by_sample[sample]
        _write_normalized(all_path, rows)
        _write_normalized(eligible_path, [row for row in rows if row[6] >= 5])
        audit = {
            "schema": "oracle-fithichip-all-interactions-v1",
            "sample": sample,
            "source_kind": "fithichip_all_interactions",
            "source_relative": f"native/{sample}.interactions_FitHiC.bed",
            "fithichip_q_filter": None,
            "merge_nearby": False,
            "eligible_min_count": 5,
            "bin_size": 5_000,
            "lower_distance": 10_000,
            "upper_distance": 50_000,
            "input_rows": len(rows),
            "retained_all_rows": len(rows),
            "retained_abundance_eligible_rows": sum(row[6] >= 5 for row in rows),
        }
        audit_path = sample_dir / f"{sample}.interactions_FitHiC.all.audit.json"
        audit_path.write_text(json.dumps(audit))
        all_files.append(str(all_path))
        eligible_files.append(str(eligible_path))
        audit_files.append(str(audit_path))

    chromsizes = tmp_path / "chrom.sizes"
    chromsizes.write_text("chr1\t100000\n")
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("")
    bedpe = tmp_path / "diff" / "union_loops.bedpe"
    support = tmp_path / "diff" / "candidate_support.tsv"
    manifest = tmp_path / "diff" / "hypothesis_universe.json"
    universe.build_universe(
        eligible_files,
        all_files,
        audit_files,
        blacklist,
        chromsizes,
        bedpe,
        support,
        manifest,
        comparison="case_vs_control",
        expected_samples=samples,
        bin_size=5_000,
        lower_distance=10_000,
        upper_distance=50_000,
        min_count=5,
        min_samples=2,
        chunk_size=1,
    )
    return samples, (loop1, loop2), by_sample, all_files, audit_files, bedpe, support, manifest


def _write_count_tables(
    tmp_path: Path, samples: list[str], loops: tuple[tuple, ...], by_sample: dict
) -> list[str]:
    outputs = []
    for sample in samples:
        source = {tuple(row[:6]): row[6] for row in by_sample[sample]}
        rows = [(*loop, source.get(tuple(loop), 0), sample) for loop in loops]
        path = tmp_path / "counts" / f"{sample}.counts.tsv"
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows, columns=COORDS + ["count", "sample"]).to_csv(
            path, sep="\t", index=False
        )
        outputs.append(str(path))
    return outputs


def test_exact_universe_keeps_two_of_four_support_and_all_source_counts(tmp_path):
    (
        samples, loops, by_sample, _all_files, audits, bedpe, support, manifest,
    ) = _universe_fixture(tmp_path)
    universe_rows = diff._load_universe(bedpe)
    assert [tuple(row) for row in universe_rows[COORDS].itertuples(index=False, name=None)] == list(loops)

    support_rows = pd.read_csv(support, sep="\t")
    assert support_rows["sample_support"].tolist() == [2, 2]
    first_counts = json.loads(support_rows.loc[0, "source_counts"])
    assert first_counts == {"s1": 8, "s2": 7, "s3": 3, "s4": 0}
    assert support_rows.loc[0, "support_samples"] == "s1,s2"
    assert (support_rows["fdr"] if "fdr" in support_rows else pd.Series()).empty

    count_files = _write_count_tables(tmp_path, samples, loops, by_sample)
    matrix, order, count_coords = diff._load_count_table(
        count_files, universe_rows, samples
    )
    assert order == samples
    contract = diff.validate_hypothesis_contract(
        bedpe,
        support,
        manifest,
        audits,
        samples,
        5,
        2,
        matrix,
        count_coords,
    )
    assert contract["condition_labels_used_for_selection"] is False
    assert contract["tolerance_bins"] == 0


def test_contract_rejects_missing_fractional_and_source_mismatch(tmp_path):
    samples, loops, by_sample, _all, audits, bedpe, support, manifest = _universe_fixture(tmp_path)
    universe_rows = diff._load_universe(bedpe)
    count_files = _write_count_tables(tmp_path, samples, loops, by_sample)

    broken = pd.read_csv(count_files[0], sep="\t").iloc[:1]
    broken.to_csv(count_files[0], sep="\t", index=False)
    with pytest.raises(ValueError, match="exact ordered"):
        diff._load_count_table(count_files, universe_rows, samples)

    count_files = _write_count_tables(tmp_path, samples, loops, by_sample)
    fractional = pd.read_csv(count_files[0], sep="\t").astype({"count": float})
    fractional.loc[0, "count"] = 1.5
    fractional.to_csv(count_files[0], sep="\t", index=False)
    with pytest.raises(ValueError, match="non-negative integers"):
        diff._load_count_table(count_files, universe_rows, samples)

    count_files = _write_count_tables(tmp_path, samples, loops, by_sample)
    missing_sample = pd.read_csv(count_files[0], sep="\t")
    missing_sample.loc[0, "sample"] = np.nan
    missing_sample.to_csv(count_files[0], sep="\t", index=False)
    with pytest.raises(ValueError, match="missing sample IDs"):
        diff._load_count_table(count_files, universe_rows, samples)

    count_files = _write_count_tables(tmp_path, samples, loops, by_sample)
    matrix, _order, count_coords = diff._load_count_table(
        count_files, universe_rows, samples
    )
    matrix.loc[matrix.index[0], "s3"] = 4
    with pytest.raises(ValueError, match="disagrees with cooler counts"):
        diff.validate_hypothesis_contract(
            bedpe, support, manifest, audits, samples, 5, 2, matrix, count_coords
        )


def test_contract_rejects_q_filtered_manifest(tmp_path):
    samples, loops, by_sample, _all, audits, bedpe, support, manifest = _universe_fixture(tmp_path)
    universe_rows = diff._load_universe(bedpe)
    count_files = _write_count_tables(tmp_path, samples, loops, by_sample)
    matrix, _order, count_coords = diff._load_count_table(
        count_files, universe_rows, samples
    )
    payload = json.loads(manifest.read_text())
    payload["fithichip_q_filter"] = 0.01
    manifest.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="fithichip_q_filter"):
        diff.validate_hypothesis_contract(
            bedpe, support, manifest, audits, samples, 5, 2, matrix, count_coords
        )


def test_universe_rejects_eligible_derivative_that_disagrees_with_all_table(tmp_path):
    samples, _loops, _by_sample, all_files, audits, _bedpe, _support, _manifest = (
        _universe_fixture(tmp_path)
    )
    eligible = tmp_path / "s1/s1.interactions_FitHiC.eligible.tsv.gz"
    frame = pd.read_csv(eligible, sep="\t")
    frame.loc[0, "score"] = 9
    frame.to_csv(eligible, sep="\t", index=False, compression="gzip")
    with pytest.raises(ValueError, match="disagree"):
        universe.build_universe(
            [
                str(tmp_path / sample / f"{sample}.interactions_FitHiC.eligible.tsv.gz")
                for sample in samples
            ],
            all_files,
            audits,
            tmp_path / "blacklist.bed",
            tmp_path / "chrom.sizes",
            tmp_path / "bad/union.bedpe",
            tmp_path / "bad/support.tsv",
            tmp_path / "bad/manifest.json",
            comparison="bad_derivative",
            expected_samples=samples,
            bin_size=5_000,
            lower_distance=10_000,
            upper_distance=50_000,
            min_count=5,
            min_samples=2,
            chunk_size=1,
        )


def test_pilot_status_and_paired_effect_schema():
    coords = pd.DataFrame([
        ("a", "chr1", 0, 5_000, "chr1", 20_000, 25_000),
        ("b", "chr1", 5_000, 10_000, "chr1", 30_000, 35_000),
    ], columns=["loop_key", *COORDS])
    normalized = pd.DataFrame(
        {"case1": [10.0, 4.0], "ctrl1": [5.0, 4.0],
         "case2": [20.0, 2.0], "ctrl2": [10.0, 8.0]},
        index=["a", "b"],
    )
    metadata = pd.DataFrame({
        "subject_id": ["d1", "d1", "d2", "d2"],
        "condition": ["case", "control", "case", "control"],
    }, index=["case1", "ctrl1", "case2", "ctrl2"])
    pairs = diff._pairing_summary(
        metadata, "subject_id", ["case1", "case2"], ["ctrl1", "ctrl2"]
    )
    status, eligible = diff.classify_analysis_status(2, 3, False)
    assert (status, eligible) == ("PILOT_UNDERPOWERED", False)
    with pytest.raises(RuntimeError, match="below"):
        diff.classify_analysis_status(2, 3, True)
    assert diff.classify_analysis_status(3, 3, True) == ("STANDARD_INFERENCE", True)

    effects = diff.build_paired_effects(
        normalized, coords, pairs, "subject_id", status
    )
    assert len(effects) == 4
    assert set(effects["pair_id"]) == {"d1", "d2"}
    assert set(effects["subject_id"]) == {"d1", "d2"}
    assert set(effects["analysis_status"]) == {"PILOT_UNDERPOWERED"}
    first = effects[(effects["loop_key"] == "a") & (effects["pair_id"] == "d1")].iloc[0]
    assert first["paired_log2_ratio"] == pytest.approx(np.log2(11 / 6))


def _project_config():
    config = yaml.safe_load((ROOT / "config/config.yaml").read_text())
    genomes = yaml.safe_load((ROOT / "config/genome.yaml").read_text())
    samples = pd.read_csv(
        ROOT / "config/samples.tsv", sep="\t", comment="#", dtype=str,
        keep_default_na=False,
    )
    return config, genomes, samples


@pytest.mark.parametrize(
    ("path", "value", "message"),
    [
        (("differential", "hypothesis_source"), "merged_calls", "hypothesis_source"),
        (("differential", "candidate_tolerance_bins"), 1, "candidate_tolerance_bins"),
        (("fithichip", "use_p2p_background"), True, "use_p2p_background"),
        (("differential", "require_publication_ready"), True, "publication minimum"),
        (("differential", "paired_by"), None, "supports paired"),
    ],
)
def test_config_rejects_ambiguous_or_underpowered_differential_modes(path, value, message):
    config, genomes, samples = _project_config()
    broken = copy.deepcopy(config)
    broken[path[0]][path[1]] = value
    with pytest.raises(ValueError, match=message):
        validate.validate_pipeline_config(broken, genomes, samples)


def test_fithichip_paths_distinguish_all_from_reporting_calls():
    text = (ROOT / "workflow/Snakefile").read_text()
    assert 'if FITHICHIP_INT_DIR == "Peak2ALL"' in text
    assert 'f"P2PBckgr_{int(_fc.get(\'use_p2p_background\', 1))}"' in text
    assert 'FITHICHIP_ALL_RESULT_DIR = "/".join(_FITHICHIP_CORE_PARTS)' in text
    assert "Differential testing uses exact native pixels with zero" in text


@pytest.mark.parametrize(
    "row",
    [
        "chr1\tbad\t5000\tchr1\t20000\t25000\t8\t0.001\t0.01\n",
        "chr1\t0\t5000\tchr1\t20000\t25000\t8\tbad\t0.01\n",
        "chr1\t0\t5000\tchr1\t20000\t25000\t8\t0.001\t0.02\n",
    ],
)
def test_reporting_filter_fails_closed_on_corrupt_or_out_of_contract_calls(tmp_path, row):
    source = tmp_path / "reporting.bed"
    source.write_text(
        "chr1\ts1\te1\tchr2\ts2\te2\tcc\tP-Value_Bias\tQ-Value_Bias\n" + row
    )
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("")
    with pytest.raises(ValueError):
        reporting_filter.filter_loops(
            source, blacklist, min_reads=5, q_threshold=0.01
        )


def test_reporting_filter_accepts_valid_header_only_zero_call_file(tmp_path):
    source = tmp_path / "zero_calls.bed"
    source.write_text(
        "chr1\ts1\te1\tchr2\ts2\te2\tcc\tP-Value_Bias\tQ-Value_Bias\n"
    )
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("")
    kept, audit = reporting_filter.filter_loops(
        source, blacklist, min_reads=5, q_threshold=0.01
    )
    assert kept.empty
    assert audit["input"] == 0
    assert audit["retained"] == 0
