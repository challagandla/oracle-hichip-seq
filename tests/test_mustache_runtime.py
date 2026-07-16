"""Focused tests for the socket-free Mustache execution adapter."""

import importlib.util
import json
import random
import sys
import threading
from pathlib import Path
from types import SimpleNamespace

import pytest


SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


runtime = _load("mustache_runtime")
wrapper = _load("mustache_balance_aware")


def _brute_force_matches(primary, mustache, resolution, tolerance_bins):
    radius = tolerance_bins * resolution
    edges = []
    for primary_key in primary:
        for mustache_key in mustache:
            if primary_key[0] != mustache_key[0] or primary_key[2] != mustache_key[2]:
                continue
            delta1 = abs(primary_key[1] - mustache_key[1])
            delta2 = abs(primary_key[3] - mustache_key[3])
            if delta1 <= radius and delta2 <= radius:
                edges.append((delta1 + delta2, primary_key, mustache_key))
    used_primary = set()
    used_mustache = set()
    matches = []
    for _distance, primary_key, mustache_key in sorted(edges):
        if primary_key in used_primary or mustache_key in used_mustache:
            continue
        used_primary.add(primary_key)
        used_mustache.add(mustache_key)
        matches.append((primary_key, mustache_key))
    return matches


def _random_loop_keys(rng, size, resolution):
    available = [
        (chrom, left * resolution, chrom, right * resolution)
        for chrom in ("chr1", "chr2")
        for left in range(0, 30)
        for right in range(left + 1, min(left + 9, 36))
    ]
    return set(rng.sample(available, min(size, len(available))))


def test_thread_backend_shares_results_and_runs_workers_concurrently():
    barrier = threading.Barrier(3)
    with runtime.ThreadManager() as manager:
        observed = manager.list()

        def worker(value):
            barrier.wait(timeout=2)
            observed.append(value)

        processes = [
            runtime.ThreadProcess(target=worker, args=(value,))
            for value in (2, 1)
        ]
        for process in processes:
            process.start()
        barrier.wait(timeout=2)
        for process in processes:
            process.join()

    assert sorted(observed) == [1, 2]
    assert all(process.exitcode == 0 for process in processes)


def test_worker_exception_is_propagated_after_other_workers_are_drained():
    completed = threading.Event()

    def fail():
        raise ValueError("broken block")

    def finish():
        completed.set()

    with pytest.raises(runtime.MustacheWorkerError, match="broken block"):
        with runtime.ThreadManager():
            failed = runtime.ThreadProcess(target=fail, name="bad-block")
            remaining = runtime.ThreadProcess(target=finish, name="good-block")
            failed.start()
            remaining.start()
            failed.join()

    assert completed.is_set()
    assert failed.exitcode == 1
    assert remaining.exitcode == 0


def test_runner_patches_only_mustache_concurrency_and_restores_globals():
    old_manager, old_process = object(), object()
    observed = {}
    fake = SimpleNamespace(Manager=old_manager, Process=old_process)

    def fake_main():
        observed["argv"] = list(sys.argv)
        observed["manager"] = fake.Manager
        observed["process"] = fake.Process
        with fake.Manager() as manager:
            values = manager.list()
            process = fake.Process(target=values.append, args=(7,))
            process.start()
            process.join()
            observed["values"] = list(values)

    fake.main = fake_main
    previous_argv = sys.argv
    runtime.run_mustache_threaded(["-f", "map.mcool"], module=fake)

    assert observed["argv"] == ["mustache", "-f", "map.mcool"]
    assert observed["manager"] is runtime.ThreadManager
    assert observed["process"] is runtime.ThreadProcess
    assert observed["values"] == [7]
    assert fake.Manager is old_manager
    assert fake.Process is old_process
    assert sys.argv is previous_argv


def test_output_is_validated_and_sorted_in_natural_genomic_order(tmp_path):
    output = tmp_path / "mustache.tsv"
    output.write_text(
        runtime.MUSTACHE_HEADER
        + "chr10\t20000\t30000\tchr10\t50000\t60000\t0.02\t3.2\n"
        + "chr2\t30000\t40000\tchr2\t80000\t90000\t0.03\t1.6\n"
        + "chr2\t10000\t20000\tchr2\t60000\t70000\t0.01\t1.6\n"
    )

    assert runtime.validate_and_sort_output(output, 10_000) == 3
    assert output.read_text().splitlines()[1:] == [
        "chr2\t10000\t20000\tchr2\t60000\t70000\t0.01\t1.6",
        "chr2\t30000\t40000\tchr2\t80000\t90000\t0.03\t1.6",
        "chr10\t20000\t30000\tchr10\t50000\t60000\t0.02\t3.2",
    ]


def test_grid_matching_exactly_matches_brute_force_greedy_randomized():
    rng = random.Random(20260716)
    resolution = 10_000
    for _iteration in range(100):
        primary = _random_loop_keys(rng, rng.randrange(0, 35), resolution)
        mustache = _random_loop_keys(rng, rng.randrange(0, 35), resolution)
        tolerance = rng.randrange(0, 4)
        observed, diagnostics = wrapper.reciprocal_anchor_matches(
            primary, mustache, resolution, tolerance
        )
        expected = _brute_force_matches(
            primary, mustache, resolution, tolerance
        )
        assert observed == expected
        assert diagnostics["grid_lookups"] == (
            len(primary) * (2 * tolerance + 1) ** 2
        )


def test_grid_matching_large_scale_guard_is_linear_in_primary_loop_count():
    resolution = 10_000
    n_loops = 20_000
    primary = {
        ("chr1", i * 3 * resolution, "chr1", (i * 3 + 100_000) * resolution)
        for i in range(n_loops)
    }
    mustache = set(primary)

    matches, diagnostics = wrapper.reciprocal_anchor_matches(
        primary, mustache, resolution, tolerance_bins=1
    )

    assert len(matches) == n_loops
    assert diagnostics["grid_lookups"] == n_loops * 9
    assert diagnostics["candidate_edges"] == n_loops
    assert diagnostics["grid_lookups"] < (len(primary) * len(mustache)) // 100


@pytest.mark.parametrize(
    "row, message",
    [
        ("chr1\t0\t9000\tchr1\t20000\t30000\t0.01\t1.6", "off the 10000-bp grid"),
        ("chr1\t0\t10000\tchr2\t20000\t30000\t0.01\t1.6", "not a valid cis loop"),
        ("chr1\t0\t10000\tchr1\t20000\t30000\t1.5\t1.6", "invalid FDR"),
        ("chr1\t0\t10000\tchr1\t20000\t30000\t0.01\tnan", "invalid detection scale"),
    ],
)
def test_invalid_caller_rows_fail_closed(tmp_path, row, message):
    output = tmp_path / "invalid.tsv"
    output.write_text(runtime.MUSTACHE_HEADER + row + "\n")

    with pytest.raises(RuntimeError, match=message):
        runtime.validate_and_sort_output(output, 10_000)


def test_wrapper_publishes_only_validated_output_and_records_backend(
    monkeypatch, tmp_path
):
    balance = tmp_path / "balance.json"
    balance.write_text(json.dumps({
        "schema": "oracle-hichip-balance-qc-v1",
        "status": "PASS",
        "weight_name": "weight",
        "resolutions": {
            "10000": {"status": "PASS", "converged": True},
        },
    }))
    primary = tmp_path / "primary.bedpe"
    primary.write_text("")
    output = tmp_path / "calls.tsv"
    status = tmp_path / "calls.status.json"

    def fake_caller(arguments):
        temporary = Path(arguments[arguments.index("-o") + 1])
        temporary.write_text(
            runtime.MUSTACHE_HEADER
            + "chr2\t20000\t30000\tchr2\t50000\t60000\t0.02\t3.2\n"
            + "chr1\t10000\t20000\tchr1\t40000\t50000\t0.01\t1.6\n"
        )

    monkeypatch.setattr(wrapper, "run_mustache_threaded", fake_caller)
    snakemake = SimpleNamespace(
        input=SimpleNamespace(
            mcool=str(tmp_path / "map.mcool"),
            balance=str(balance),
            primary=str(primary),
        ),
        output=SimpleNamespace(tsv=str(output), status=str(status)),
        params=SimpleNamespace(res=10_000, comparison_tolerance_bins=1),
        threads=2,
        wildcards=SimpleNamespace(sample="sample_1"),
        log=[str(tmp_path / "mustache.log")],
    )

    wrapper.main(snakemake)

    payload = json.loads(status.read_text())
    assert payload["status"] == "PASS"
    assert payload["execution_backend"] == runtime.EXECUTION_BACKEND
    assert payload["n_output_rows"] == 2
    assert output.read_text().splitlines()[1].startswith("chr1\t")
    assert not list(tmp_path.glob(".*.mustache-*.tmp"))


def test_wrapper_does_not_publish_partial_output_after_caller_error(
    monkeypatch, tmp_path
):
    balance = tmp_path / "balance.json"
    balance.write_text(json.dumps({
        "status": "PASS",
        "weight_name": "weight",
        "resolutions": {
            "10000": {"status": "PASS", "converged": True},
        },
    }))
    primary = tmp_path / "primary.bedpe"
    primary.write_text("")
    output = tmp_path / "calls.tsv"
    status = tmp_path / "calls.status.json"

    def broken_caller(arguments):
        Path(arguments[arguments.index("-o") + 1]).write_text(
            runtime.MUSTACHE_HEADER
        )
        raise RuntimeError("caller failed")

    monkeypatch.setattr(wrapper, "run_mustache_threaded", broken_caller)
    snakemake = SimpleNamespace(
        input=SimpleNamespace(
            mcool=str(tmp_path / "map.mcool"),
            balance=str(balance),
            primary=str(primary),
        ),
        output=SimpleNamespace(tsv=str(output), status=str(status)),
        params=SimpleNamespace(res=10_000, comparison_tolerance_bins=1),
        threads=2,
        wildcards=SimpleNamespace(sample="sample_1"),
        log=[str(tmp_path / "mustache.log")],
    )

    with pytest.raises(RuntimeError, match="caller failed"):
        wrapper.main(snakemake)

    assert not output.exists()
    assert not status.exists()
    assert not list(tmp_path.glob(".*.mustache-*.tmp"))
