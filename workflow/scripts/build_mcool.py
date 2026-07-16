"""Build, audit and safely publish a balanced multi-resolution cooler."""
import gzip
import os
import shlex
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from balance_qc import write_balance_tsv  # noqa: E402
from balance_utils import inspect_mcool_balance, sanitize_nonpass_weights  # noqa: E402
from utils import setup_logging, write_json  # noqa: E402


def build_zoomify_command(
    input_cool: str,
    output_mcool: str,
    resolutions_bp: list[int],
    threads: int,
    *,
    weight_name: str,
    ignore_diags: int,
    min_nnz: int,
    mad_max: int,
    tolerance: float,
    max_iterations: int,
    blacklist: str | None = None,
) -> list[str]:
    balance_args_parts = [
        "--nproc", str(threads),
        "--name", str(weight_name),
        "--mad-max", str(mad_max),
        "--min-nnz", str(min_nnz),
        "--ignore-diags", str(ignore_diags),
        "--tol", str(tolerance),
        "--max-iters", str(max_iterations),
        "--convergence-policy", "store_final",
    ]
    if blacklist:
        balance_args_parts.extend(["--blacklist", blacklist])
    balance_args = shlex.join(balance_args_parts)
    return [
        "cooler", "zoomify",
        "--nproc", str(threads),
        "--balance",
        "--balance-args", balance_args,
        "--resolutions", ",".join(str(value) for value in resolutions_bp),
        "-o", output_mcool,
        input_cool,
    ]


def materialize_plaintext_blacklist(
    source: str | Path,
    destination: str | Path,
) -> Path:
    """Copy a BED blacklist to plain UTF-8 text for ``cooler balance``.

    Reference blacklists are commonly distributed as gzip-compressed BED files,
    while Cooler requires an uncompressed three-column BED for ``--blacklist``.
    Detect compression by file signature rather than filename so ``.bgz`` and
    extensionless inputs are handled correctly as well.
    """
    source_path = Path(source)
    destination_path = Path(destination)
    with source_path.open("rb") as handle:
        is_gzip = handle.read(2) == b"\x1f\x8b"

    opener = gzip.open if is_gzip else source_path.open
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    if is_gzip:
        source_handle = opener(source_path, "rt", encoding="utf-8")
    else:
        source_handle = opener("r", encoding="utf-8")
    with source_handle as src, destination_path.open(
        "w", encoding="utf-8", newline=""
    ) as dst:
        shutil.copyfileobj(src, dst)
    return destination_path


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    output_mcool = Path(snakemake.output.mcool)
    output_mcool.parent.mkdir(parents=True, exist_ok=True)
    attempted = output_mcool.with_name(
        f".{output_mcool.name}.{uuid.uuid4().hex}.attempted"
    )
    balance_blacklist = attempted.with_name(f"{attempted.name}.blacklist.bed")
    resolutions = sorted({int(value) for value in snakemake.params.resolutions_bp})
    requested = {
        "weight_name": str(snakemake.params.weight_name),
        "ignore_diags": int(snakemake.params.ignore_diags),
        "min_nnz": int(snakemake.params.min_nnz),
        "mad_max": int(snakemake.params.mad_max),
        "tolerance": float(snakemake.params.tolerance),
        "max_iterations": int(snakemake.params.max_iterations),
        "convergence_policy": "store_final_then_remove_nonpass_from_published_mcool",
    }
    try:
        materialize_plaintext_blacklist(
            snakemake.input.blacklist, balance_blacklist
        )
        command = build_zoomify_command(
            str(snakemake.input.cool), str(attempted), resolutions,
            int(getattr(snakemake, "threads", 1)),
            weight_name=requested["weight_name"],
            ignore_diags=requested["ignore_diags"],
            min_nnz=requested["min_nnz"],
            mad_max=requested["mad_max"],
            tolerance=requested["tolerance"],
            max_iterations=requested["max_iterations"],
            blacklist=str(balance_blacklist),
        )
        with Path(snakemake.log[0]).open("a") as log:
            subprocess.run(command, check=True, stdout=log, stderr=log)
        report = inspect_mcool_balance(
            attempted, resolutions, requested["weight_name"]
        )
        report["sample"] = snakemake.wildcards.sample
        report["mcool"] = str(output_mcool)
        report["requested_balance_parameters"] = requested
        sanitize_nonpass_weights(attempted, report)
        os.replace(attempted, output_mcool)
        write_json(report, snakemake.output.json)
        write_balance_tsv(report, snakemake.output.tsv)
    finally:
        attempted.unlink(missing_ok=True)
        balance_blacklist.unlink(missing_ok=True)


if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
