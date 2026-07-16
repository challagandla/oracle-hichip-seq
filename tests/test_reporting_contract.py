"""Regression tests for reporting, provenance, and transparent text inputs."""
import copy
import gzip
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
        f"reporting_contract_{name}", SCRIPTS / f"{name}.py"
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


utils = _load("utils")
filter_blacklist = _load("filter_cool_blacklist")
consensus = _load("build_consensus_loops")
apa = _load("apa_plot")
restriction = _load("restriction_fragment_qc")
figures = _load("figures")
multiqc = _load("multiqc_content")
provenance = _load("provenance_manifest")
validate_config = _load("validate_config")


@pytest.mark.parametrize("compressed", [False, True])
def test_open_text_auto_uses_magic_bytes_not_suffix(tmp_path, compressed):
    path = tmp_path / ("actually_plain.gz" if not compressed else "compressed.bed")
    content = "chr1\t10\t20\n"
    if compressed:
        with gzip.open(path, "wt") as handle:
            handle.write(content)
    else:
        path.write_text(content)

    with utils.open_text_auto(path) as handle:
        assert handle.read() == content

    assert filter_blacklist._read_blacklist(path)["chr1"][0].tolist() == [10]
    assert consensus._read_blacklist(path)["chr1"][0].tolist() == [10]
    assert apa._read_interval_index(path)["chr1"][0].tolist() == [10]


def test_loop_and_restriction_readers_accept_misleading_suffixes(tmp_path):
    loops = tmp_path / "loops.bedpe"
    with gzip.open(loops, "wt") as handle:
        handle.write("chr1\t0\t5000\tchr1\t10000\t15000\n")
    assert len(utils.load_loops_bedpe(loops)) == 1

    pairs = tmp_path / "contacts.pairs.gz"
    pairs.write_text(
        "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type rfrag1 rfrag2\n"
        "r1\tchr1\t100\tchr1\t500\t+\t-\tUU\t1\t2\n"
    )
    report = restriction.summarise_pairs(pairs)
    assert report["population"] == "post_dedup_pre_contact_filter_UU_pairs"
    assert report["total_deduplicated_uu_pairs"] == 1


def test_required_figure_inputs_fail_closed(tmp_path):
    malformed = tmp_path / "broken.json"
    malformed.write_text("not-json")
    with pytest.raises(ValueError, match="not valid JSON"):
        figures._required_json(malformed, "test report", {"status"})

    incomplete = tmp_path / "incomplete.tsv"
    incomplete.write_text("other\n1\n")
    with pytest.raises(ValueError, match="lacks columns"):
        figures._required_tsv(incomplete, "test table", {"required"})

    malformed_loops = tmp_path / "loops.tsv"
    malformed_loops.write_text(
        "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\n"
        "chr1\tnot-a-number\t5000\tchr1\t10000\t15000\n"
    )
    with pytest.raises(ValueError, match="nonnumeric start1"):
        figures._required_pipeline_loops(malformed_loops, "pipeline loops")


def test_apa_status_and_matrix_must_agree():
    zero = np.zeros((5, 5))
    nonzero = zero.copy()
    nonzero[2, 2] = 1
    assert figures.validate_apa_matrix_contract(
        {"status": "NOT_ASSESSED"}, zero, "sample"
    ) == "NOT_ASSESSED"
    assert figures.validate_apa_matrix_contract(
        {"status": "DESCRIPTIVE"}, nonzero, "sample"
    ) == "DESCRIPTIVE"
    with pytest.raises(ValueError, match="NOT_ASSESSED"):
        figures.validate_apa_matrix_contract(
            {"status": "NOT_ASSESSED"}, nonzero, "sample"
        )


def test_figures_use_configured_gates_and_primary_reporting_role():
    assert figures.qc_gate_values({
        "valid_pair_yield_pct": 31,
        "duplicate_pct_max": 44,
        "cis_fraction_min": 0.73,
    }) == pytest.approx((31, 44, 73))
    with pytest.raises(ValueError, match="missing"):
        figures.qc_gate_values({"duplicate_pct_max": 44})

    values = pd.Series([20.0, 50.0, 80.0, np.nan])
    assert figures.qc_gate_failures(values, minimum=50).tolist() == [
        True, False, False, False
    ]
    assert figures.qc_gate_failures(values, maximum=50).tolist() == [
        False, False, True, False
    ]

    library = pd.DataFrame(
        {"report_role": ["primary", "demonstration"]}, index=["real", "demo"]
    )
    assert figures.primary_reporting_libraries(library).index.tolist() == ["real"]


def _write_differential_bundle(root: Path, status: str = "PILOT_UNDERPOWERED"):
    base = root / "diff" / "case_vs_control"
    base.mkdir(parents=True)
    pd.DataFrame([{
        "loop_key": "abc123", "chrom1": "chr1", "start1": 0, "end1": 5000,
        "chrom2": "chr1", "start2": 10000, "end2": 15000,
        "padj": 0.04, "log2FoldChange": 1.2, "lfcSE": 0.4,
        "analysis_status": status,
    }]).to_csv(base / "differential_loops.tsv", sep="\t", index=False)
    design = {
        "comparison": "case_vs_control",
        "analysis_status": status,
        "n_complete_pairs": 2 if status == "PILOT_UNDERPOWERED" else 3,
        "publication_eligible": status == "STANDARD_INFERENCE",
        "publication_min_complete_pairs": 3,
        "paired_subjects": ["D1", "D2"],
        "candidate_loops": 1,
        "tested_loops": 1,
    }
    (base / "design.json").write_text(json.dumps(design))
    pd.DataFrame([
        {
            "loop_key": "abc123", "pairing_factor": "subject_id", "pair_id": pair,
            "chrom1": "chr1", "start1": 0, "end1": 5000,
            "chrom2": "chr1", "start2": 10000, "end2": 15000,
            "subject_id": pair, "case_sample": f"case{i + 1}",
            "control_sample": f"ctrl{i + 1}", "case_normalized_count": 20.0 + i,
            "control_normalized_count": 10.0 + i,
            "paired_log2_ratio": 0.93 - i / 10, "analysis_status": status,
        }
        for i, pair in enumerate(design["paired_subjects"])
    ]).to_csv(base / "paired_effects.tsv", sep="\t", index=False)
    return base


def test_differential_figure_contract_requires_pilot_and_paired_effects(tmp_path):
    base = _write_differential_bundle(tmp_path)
    result, design, paired = figures.load_differential_bundle(
        "case_vs_control", tmp_path
    )
    assert design["analysis_status"] == "PILOT_UNDERPOWERED"
    assert result["lfcSE"].iloc[0] == pytest.approx(0.4)
    assert paired["paired_log2_ratio"].iloc[0] == pytest.approx(0.93)

    broken = pd.read_csv(base / "paired_effects.tsv", sep="\t")
    broken["analysis_status"] = "STANDARD_INFERENCE"
    broken.to_csv(base / "paired_effects.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match="disagrees"):
        figures.load_differential_bundle("case_vs_control", tmp_path)


def test_differential_figure_contract_rejects_truncation_and_bad_bool(tmp_path):
    truncated_root = tmp_path / "truncated"
    base = _write_differential_bundle(truncated_root)
    paired = pd.read_csv(base / "paired_effects.tsv", sep="\t").iloc[:1]
    paired.to_csv(base / "paired_effects.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match="paired effects has 1 rows"):
        figures.load_differential_bundle("case_vs_control", truncated_root)

    duplicate_root = tmp_path / "duplicate"
    base = _write_differential_bundle(duplicate_root)
    paired = pd.read_csv(base / "paired_effects.tsv", sep="\t")
    pd.concat([paired, paired.iloc[[0]]], ignore_index=True).to_csv(
        base / "paired_effects.tsv", sep="\t", index=False
    )
    with pytest.raises(ValueError, match="duplicate"):
        figures.load_differential_bundle("case_vs_control", duplicate_root)

    bool_root = tmp_path / "bad_bool"
    base = _write_differential_bundle(bool_root)
    design = json.loads((base / "design.json").read_text())
    design["publication_eligible"] = "false"
    (base / "design.json").write_text(json.dumps(design))
    with pytest.raises(ValueError, match="JSON boolean"):
        figures.load_differential_bundle("case_vs_control", bool_root)


@pytest.mark.parametrize(
    ("column", "bad_value", "message"),
    [
        ("padj", -0.01, "padj below"),
        ("padj", 1.01, "padj above"),
        ("padj", np.inf, "non-finite padj"),
        ("log2FoldChange", np.inf, "non-finite log2FoldChange"),
        ("lfcSE", -0.1, "lfcSE below"),
        ("lfcSE", np.inf, "non-finite lfcSE"),
    ],
)
def test_differential_figure_contract_rejects_invalid_model_values(
    tmp_path, column, bad_value, message
):
    root = tmp_path / column / str(bad_value)
    base = _write_differential_bundle(root)
    result = pd.read_csv(base / "differential_loops.tsv", sep="\t")
    result.loc[0, column] = bad_value
    result.to_csv(base / "differential_loops.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match=message):
        figures.load_differential_bundle("case_vs_control", root)


def test_stripe_table_requires_caller_schema_and_positive_length(tmp_path):
    wrong = tmp_path / "wrong.tsv"
    wrong.write_text("length\n10000\n")
    with pytest.raises(ValueError, match="lacks columns"):
        figures._required_stripe_table(wrong, "stripe output")

    columns = list(figures.STRIPE_COLUMNS)
    invalid = tmp_path / "invalid.tsv"
    pd.DataFrame([{column: 1 for column in columns}]).assign(length=-1).to_csv(
        invalid, sep="\t", index=False
    )
    with pytest.raises(ValueError, match="non-positive"):
        figures._required_stripe_table(invalid, "stripe output")


def test_multiqc_exposes_both_frips_and_differential_status():
    report = {
        "sample": "sample_1", "valid_pair_yield_pct": 40.0,
        "duplicate_pct": 20.0, "cis_fraction": 0.8,
        "sample_peak_frip": 0.5, "consensus_peak_frip": 0.4,
        "n_sample_peaks": 100, "n_consensus_peaks": 80, "n_loops": 200,
        "apa_matched_effect": 2.1, "apa_corner_ratio_near_distance": 1.8,
        "apa_status": "DESCRIPTIVE", "apa_normalization": "ICE-balanced",
        "overall_status": "PASS", "thresholds": {"hicrep_scc_min": 0.7},
    }
    payload = multiqc.loop_qc_payload(report)
    row = payload["data"]["sample_1"]
    assert row["sample_peak_frip"] == pytest.approx(0.5)
    assert row["consensus_peak_frip"] == pytest.approx(0.4)
    assert "matched" in payload["headers"]["apa_matched_effect"]["title"].lower()
    assert "+/-2 bins" in payload["headers"]["apa_corner_ratio"]["description"]

    differential = multiqc.differential_payload({
        "comparison": "case_vs_control", "analysis_status": "PILOT_UNDERPOWERED",
        "n_complete_pairs": 2, "publication_eligible": False,
        "publication_min_complete_pairs": 3,
    })
    assert differential["data"]["case_vs_control"]["publication_eligible"] is False
    assert "exploratory" in differential["description"]


def test_provenance_records_only_portable_package_fields(tmp_path):
    records = provenance.normalize_conda_records([{
        "name": "pandas", "version": "2.2.3", "build_string": "py311_0",
        "channel": "conda-forge", "platform": "linux-64",
        "prefix": "/host/specific/path", "base_url": "https://example.invalid",
        "build_number": 0,
    }])
    assert set(records[0]) == set(provenance.PORTABLE_PACKAGE_FIELDS)
    assert "/host/specific/path" not in json.dumps(records)

    spec = tmp_path / "pandas.yaml"
    spec.write_text("name: test\ndependencies:\n  - python=3.11\n")
    cache = tmp_path / "conda"
    cache.mkdir()
    captured = cache / "abc_.yaml"
    captured.write_bytes(spec.read_bytes())
    prefix = captured.with_suffix("")
    (prefix / "conda-meta").mkdir(parents=True)
    assert provenance.matching_snakemake_conda_prefix(spec, cache) == prefix


def test_provenance_hashes_named_reference_assets(tmp_path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGT\n")
    records = provenance.reference_records({"fasta": str(fasta)})
    assert records["fasta"]["sha256"] == provenance.sha256(fasta)
    assert records["fasta"]["size_bytes"] == fasta.stat().st_size


def test_rules_declare_reporting_and_reference_contracts():
    peaks = (ROOT / "workflow/rules/04_peaks.smk").read_text()
    figures_rule = (ROOT / "workflow/rules/11_figures.smk").read_text()
    provenance_rule = (ROOT / "workflow/rules/12_provenance.smk").read_text()
    assert "-nonamecheck" in peaks
    assert "primary_autosomes_chrX" in peaks
    assert "sample_peak_frip" in peaks and "consensus_peak_frip" in peaks
    assert "qc_thresholds = config[\"qc_thresholds\"]" in figures_rule
    assert "paired_effects.tsv" in figures_rule
    assert "reference_alignment_indexes = BWA_INDEX_FILES" in provenance_rule
    assert "reference_contract = PROVENANCE_REFERENCE_CONTRACT" in provenance_rule
    assert "hypothesis_universes" in provenance_rule
    assert "interactions_FitHiC.all.audit.json" in provenance_rule


def test_apa_reconciliation_is_intentionally_distinct_from_differential_grid():
    config = yaml.safe_load((ROOT / "config/config.yaml").read_text())
    genomes = yaml.safe_load((ROOT / "config/genome.yaml").read_text())
    samples = pd.read_csv(
        ROOT / "config/samples.tsv", sep="\t", comment="#", dtype=str,
        keep_default_na=False,
    )
    assert config["apa"]["candidate_tolerance_bins"] == 1
    assert config["differential"]["candidate_tolerance_bins"] == 0
    rule = (ROOT / "workflow/rules/06_loop_qc.smk").read_text()
    assert 'tolerance_bins = config["apa"]["candidate_tolerance_bins"]' in rule
    assert 'candidate_grid_bin_size_bp = config["fithichip"]["bin_size"]' in rule

    broken = copy.deepcopy(config)
    broken["apa"]["candidate_tolerance_bins"] = 0
    with pytest.raises(ValueError, match="apa.candidate_tolerance_bins"):
        validate_config.validate_pipeline_config(broken, genomes, samples)
