"""Run expected-cis or insulation with balance-aware normalization fallback."""
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from balance_utils import load_balance_report, resolution_balance  # noqa: E402
from utils import setup_logging  # noqa: E402


def build_command(
    kind: str,
    matrix: str,
    output: str,
    view: str,
    resolution: int,
    threads: int,
    ignore_diags: int,
    use_balanced: bool,
    weight_name: str = "weight",
    window: int | None = None,
) -> list[str]:
    if kind == "expected_cis":
        cmd = [
            "cooltools", "expected-cis", "-p", str(threads),
            "--view", view, "--ignore-diags", str(ignore_diags),
        ]
    elif kind == "insulation":
        if window is None:
            raise ValueError("insulation requires a window size")
        cmd = [
            "cooltools", "insulation", "-p", str(threads),
            "--view", view, "--ignore-diags", str(ignore_diags),
        ]
    else:
        raise ValueError(f"Unknown cooltools matrix-QC kind {kind!r}")

    # Pass an actual empty subprocess argument for raw mode. Shell quoting such
    # as --clr-weight-name '' is not equivalent when commands are constructed
    # as a list; the empty string itself is the cooltools raw-count sentinel.
    cmd.extend(["--clr-weight-name", weight_name if use_balanced else ""])
    cmd.extend([f"{matrix}::resolutions/{resolution}"])
    if kind == "insulation":
        cmd.append(str(window))
    cmd.extend(["-o", output])
    return cmd


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    kind = str(snakemake.params.kind)
    resolution = int(snakemake.params.res)
    decision = resolution_balance(
        load_balance_report(snakemake.input.balance), resolution
    )
    output = Path(snakemake.output.tsv)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".cooltools.tmp")
    command = build_command(
        kind=kind,
        matrix=str(snakemake.input.mcool),
        output=str(temporary),
        view=str(snakemake.input.view),
        resolution=resolution,
        threads=int(getattr(snakemake, "threads", 1)),
        ignore_diags=int(snakemake.params.ignore_diags),
        use_balanced=bool(decision["use_balanced"]),
        weight_name=str(decision["weight_name"]),
        window=(int(snakemake.params.window) if kind == "insulation" else None),
    )
    with open(snakemake.log[0], "a") as log:
        log.write(
            f"normalization={decision['normalization']} "
            f"balance_status={decision['status']} resolution_bp={resolution}\n"
        )
        subprocess.run(command, check=True, stdout=log, stderr=log)

    table = pd.read_csv(temporary, sep="\t")
    table["normalization"] = decision["normalization"]
    table["balance_status"] = decision["status"]
    table["balance_converged"] = decision.get("converged")
    table["resolution_bp"] = resolution
    table.to_csv(output, sep="\t", index=False)
    temporary.unlink(missing_ok=True)


if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
