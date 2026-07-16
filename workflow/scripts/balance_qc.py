"""Write stable JSON and TSV QC from cooler balancing dataset attributes."""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from balance_utils import inspect_mcool_balance  # noqa: E402
from utils import setup_logging, write_json  # noqa: E402


TSV_COLUMNS = [
    "resolution_bp", "status", "weight_present", "weight_published",
    "converged", "variance",
    "tolerance", "cis_only", "divisive_weights", "ignore_diags", "mad_max",
    "min_count", "min_nnz", "scale",
]


def write_balance_tsv(report: dict, path: str | Path) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t")
        writer.writeheader()
        for key in report["configured_resolutions_bp"]:
            entry = report["resolutions"][str(key)]
            params = entry.get("parameters", {})
            writer.writerow({
                "resolution_bp": key,
                "status": entry["status"],
                "weight_present": entry["weight_present"],
                "weight_published": entry.get("weight_published", False),
                "converged": entry["converged"],
                "variance": entry["variance"],
                "tolerance": entry["tolerance"],
                **{name: params.get(name) for name in TSV_COLUMNS[7:]},
            })


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    report = inspect_mcool_balance(
        snakemake.input.mcool,
        list(snakemake.params.resolutions_bp),
        str(snakemake.params.weight_name),
    )
    report["sample"] = snakemake.wildcards.sample
    write_json(report, snakemake.output.json)
    write_balance_tsv(report, snakemake.output.tsv)


if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
