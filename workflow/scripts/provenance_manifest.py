"""Write a portable run manifest without copying raw data or environment paths."""
import hashlib
import json
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


PORTABLE_PACKAGE_FIELDS = ("name", "version", "build_string", "channel", "platform")


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json_sha256(value) -> str:  # type: ignore[no-untyped-def]
    """Hash the effective configuration independently of YAML formatting."""
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), default=str
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def git_value(repository: str, *args: str) -> str | None:
    proc = subprocess.run(
        ["git", "-C", repository, *args], capture_output=True, text=True, check=False
    )
    return proc.stdout.strip() if proc.returncode == 0 else None


def normalize_conda_records(records: list[dict]) -> list[dict[str, str | None]]:
    """Strip host paths and solver metadata from ``conda list --json`` output."""
    portable = []
    for record in records:
        if not isinstance(record, dict) or not record.get("name"):
            raise ValueError("conda package records must be objects with a name")
        portable.append({field: record.get(field) for field in PORTABLE_PACKAGE_FIELDS})
    return sorted(
        portable,
        key=lambda item: (
            str(item["name"]), str(item["version"]), str(item["build_string"])
        ),
    )


def conda_package_records(
    conda_executable: str, *, prefix: str | Path | None = None,
    environment_name: str | None = None,
) -> list[dict[str, str | None]]:
    """Return portable resolved package builds for one installed environment."""
    if (prefix is None) == (environment_name is None):
        raise ValueError("provide exactly one of prefix or environment_name")
    command = [conda_executable, "list", "--json"]
    if prefix is not None:
        command += ["--prefix", str(prefix)]
    else:
        command += ["--name", str(environment_name)]
    proc = subprocess.run(command, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        target = str(prefix if prefix is not None else environment_name)
        raise RuntimeError(
            f"could not resolve installed packages for {target}: {proc.stderr.strip()}"
        )
    try:
        records = json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError("conda list returned invalid JSON") from exc
    if not isinstance(records, list) or not records:
        raise RuntimeError("conda list returned no installed package records")
    return normalize_conda_records(records)


def matching_snakemake_conda_prefix(
    environment_yaml: str | Path, cache_directory: str | Path
) -> Path:
    """Find the installed Snakemake prefix whose captured YAML is byte-identical."""
    specification = Path(environment_yaml)
    expected_hash = sha256(specification)
    matches = []
    for captured in Path(cache_directory).glob("*_.yaml"):
        if sha256(captured) != expected_hash:
            continue
        prefix = captured.with_suffix("")
        if (prefix / "conda-meta").is_dir():
            matches.append(prefix)
    if len(matches) != 1:
        raise RuntimeError(
            f"expected one installed Snakemake environment for {specification}; "
            f"found {len(matches)} in {cache_directory}"
        )
    return matches[0]


def reference_records(reference_contract: dict[str, str]) -> dict[str, dict]:
    """Hash every configured analysis reference under a stable semantic name."""
    records = {}
    for name, configured_path in sorted(reference_contract.items()):
        path = Path(configured_path)
        if not path.is_file():
            raise FileNotFoundError(f"required reference {name!r} is missing: {path}")
        records[str(name)] = {
            "configured_path": str(configured_path),
            "size_bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
    if not records:
        raise ValueError("reference provenance contract is empty")
    return records


def resolved_environment_records(
    environment_yamls: list[str],
    runner_yaml: str,
    cache_directory: str,
    conda_executable: str,
    runner_environment: str,
) -> dict[str, dict]:
    """Capture exact builds without embedding machine-specific prefixes."""
    records: dict[str, dict] = {
        "runner": {
            "definition": Path(runner_yaml).name,
            "definition_sha256": sha256(runner_yaml),
            "packages": conda_package_records(
                conda_executable, environment_name=runner_environment
            ),
        }
    }
    for specification in sorted(environment_yamls):
        prefix = matching_snakemake_conda_prefix(specification, cache_directory)
        key = f"workflow/envs/{Path(specification).name}"
        records[key] = {
            "definition": key,
            "definition_sha256": sha256(specification),
            "packages": conda_package_records(conda_executable, prefix=prefix),
        }
    return records


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    contract_files = [
        str(snakemake.input.samples), str(snakemake.input.parameters),
        str(snakemake.input.genome), str(snakemake.input.runner),
        *list(snakemake.input.envs),
    ]
    small_outputs = [
        str(snakemake.input.multiqc), str(snakemake.input.figures),
        *list(snakemake.input.figure_files), *list(snakemake.input.loop_qc),
        *list(snakemake.input.balance_qc), *list(snakemake.input.blacklist_qc),
        *list(snakemake.input.loop_call_audits),
        *list(snakemake.input.oracle), *list(snakemake.input.stripes),
        *list(snakemake.input.differential),
        *list(snakemake.input.hypothesis_universes),
    ]
    balance_status = {}
    for path in snakemake.input.balance_qc:
        report = json.loads(Path(path).read_text())
        balance_status[str(report.get("sample", Path(path).stem))] = {
            "status": report.get("status", "NOT_ASSESSED"),
            "nonconverged_resolutions_bp": report.get(
                "nonconverged_resolutions_bp", []
            ),
            "missing_resolutions_bp": report.get("missing_resolutions_bp", []),
        }

    references = reference_records(dict(snakemake.params.reference_contract))
    conda_executable = (
        os.environ.get("CONDA_EXE")
        or shutil.which("conda")
        or shutil.which("mamba")
    )
    if not conda_executable:
        raise RuntimeError("conda executable is unavailable; resolved builds cannot be recorded")
    environments = resolved_environment_records(
        list(snakemake.input.envs),
        str(snakemake.input.runner),
        str(snakemake.params.conda_cache),
        conda_executable,
        str(snakemake.params.runner_environment),
    )

    repository = str(snakemake.params.repository)
    status = git_value(repository, "status", "--porcelain", "--untracked-files=no")
    effective_config = json.loads(
        json.dumps(snakemake.params.effective_config, default=str)
    )
    manifest = {
        "schema": "oracle-hichip-run-manifest-v2",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "pipeline": {
            "commit": git_value(repository, "rev-parse", "HEAD"),
            "branch": git_value(repository, "branch", "--show-current"),
            "tracked_worktree_clean": status == "",
        },
        "effective_config": effective_config,
        "effective_config_sha256": canonical_json_sha256(effective_config),
        "input_contract_sha256": {path: sha256(path) for path in contract_files},
        "references": references,
        "reference_sha256": {
            name: record["sha256"] for name, record in references.items()
        },
        "resolved_environments": environments,
        "report_sha256": {path: sha256(path) for path in small_outputs},
        "balance_qc": balance_status,
        "multiqc": str(snakemake.input.multiqc),
        "note": (
            "Raw sequencing data and large derived matrices are referenced by the "
            "Snakemake DAG and are intentionally not copied or hashed here. Reference "
            "assets, environment definitions, and portable resolved package builds "
            "are recorded for reproducibility."
        ),
    }
    output = Path(snakemake.output.json)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2) + "\n")
    log = Path(snakemake.log[0])
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text(
        f"Wrote provenance v2 with {len(references)} references and "
        f"{len(environments)} resolved environments to {output}\n"
    )


if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
