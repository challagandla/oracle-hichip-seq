"""Shared cooler-balance inspection and downstream normalization decisions."""
import json
from pathlib import Path


BALANCE_STATUSES = frozenset({"PASS", "WARN", "NOT_ASSESSED"})


def _json_scalar(value):
    """Convert h5py/NumPy attribute scalars to stable JSON values."""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if hasattr(value, "tolist"):
        return value.tolist()
    return str(value)


def inspect_mcool_balance(
    mcool_path: str | Path,
    resolutions_bp: list[int],
    weight_name: str = "weight",
) -> dict:
    """Read convergence metadata from each configured ``bins/weight`` dataset.

    Dataset presence is not evidence of usable ICE weights. Cooler writes the
    dataset even when iterative correction stops without convergence, recording
    the authoritative result in the dataset's ``converged`` attribute.
    """
    import h5py

    resolution_reports: dict[str, dict] = {}
    with h5py.File(mcool_path, "r") as handle:
        for resolution in sorted({int(value) for value in resolutions_bp}):
            dataset_path = f"resolutions/{resolution}/bins/{weight_name}"
            entry = {
                "resolution_bp": resolution,
                "weight_dataset": f"/{dataset_path}",
                "weight_present": False,
                "weight_published": False,
                "status": "NOT_ASSESSED",
                "converged": None,
                "variance": None,
                "tolerance": None,
                "parameters": {},
            }
            if dataset_path not in handle:
                entry["note"] = "Balancing weight dataset is missing."
                resolution_reports[str(resolution)] = entry
                continue

            dataset = handle[dataset_path]
            attrs = {
                str(key): _json_scalar(value)
                for key, value in sorted(dataset.attrs.items())
            }
            converged = attrs.get("converged")
            if converged is not None:
                converged = bool(converged)
            entry.update({
                "weight_present": True,
                "weight_published": True,
                "converged": converged,
                "variance": attrs.get("var"),
                "tolerance": attrs.get("tol"),
                "parameters": {
                    key: value
                    for key, value in attrs.items()
                    if key not in {"converged", "var", "tol"}
                },
            })
            if converged is True:
                entry["status"] = "PASS"
            elif converged is False:
                entry["status"] = "WARN"
                entry["note"] = "Stored balancing weights did not converge."
            else:
                entry["note"] = "Balancing convergence attribute is missing."
            resolution_reports[str(resolution)] = entry

    missing = [
        int(key) for key, value in resolution_reports.items()
        if value["status"] == "NOT_ASSESSED"
    ]
    nonconverged = [
        int(key) for key, value in resolution_reports.items()
        if value["status"] == "WARN"
    ]
    converged = [
        int(key) for key, value in resolution_reports.items()
        if value["status"] == "PASS"
    ]
    if missing:
        overall_status, passed = "NOT_ASSESSED", None
    elif nonconverged:
        overall_status, passed = "WARN", None
    else:
        overall_status, passed = "PASS", True

    return {
        "schema": "oracle-hichip-balance-qc-v1",
        "mcool": str(mcool_path),
        "weight_name": weight_name,
        "configured_resolutions_bp": sorted({int(value) for value in resolutions_bp}),
        "status": overall_status,
        "pass": passed,
        "n_configured": len(resolution_reports),
        "n_converged": len(converged),
        "n_nonconverged": len(nonconverged),
        "n_missing": len(missing),
        "converged_resolutions_bp": converged,
        "nonconverged_resolutions_bp": nonconverged,
        "missing_resolutions_bp": missing,
        "resolutions": resolution_reports,
    }


def sanitize_nonpass_weights(mcool_path: str | Path, report: dict) -> list[int]:
    """Delete every configured weight dataset that is not proven converged.

    Cooler and hicmatrix clients commonly use a column named ``weight`` without
    checking its ``converged`` attribute.  The attempted-balance report retains
    those attributes, while the published mcool fails safe: only PASS weights
    remain available for implicit consumers.
    """
    import h5py

    removed: list[int] = []
    with h5py.File(mcool_path, "r+") as handle:
        for key in report.get("configured_resolutions_bp", []):
            entry = report["resolutions"][str(int(key))]
            dataset_path = str(entry.get("weight_dataset", "")).lstrip("/")
            if entry.get("status") == "PASS":
                entry["weight_published"] = dataset_path in handle
                continue
            if dataset_path and dataset_path in handle:
                del handle[dataset_path]
                removed.append(int(key))
            entry["weight_published"] = False

    report["published_weight_policy"] = (
        "Only weights with converged=true are retained in the published mcool."
    )
    report["removed_weight_resolutions_bp"] = removed
    return removed


def load_balance_report(path: str | Path) -> dict:
    report = json.loads(Path(path).read_text())
    status = str(report.get("status", "NOT_ASSESSED"))
    if status not in BALANCE_STATUSES:
        raise ValueError(f"Unknown cooler-balance status {status!r}")
    return report


def resolution_balance(report: dict, resolution_bp: int) -> dict:
    """Return one resolution's state, failing closed to raw/NOT_ASSESSED."""
    resolution = int(resolution_bp)
    entry = dict(report.get("resolutions", {}).get(str(resolution), {}))
    status = str(entry.get("status", "NOT_ASSESSED"))
    if status not in BALANCE_STATUSES:
        raise ValueError(
            f"Unknown cooler-balance status {status!r} at {resolution} bp"
        )
    use_balanced = status == "PASS"
    return {
        **entry,
        "resolution_bp": resolution,
        "status": status,
        "weight_name": str(report.get("weight_name", "weight")),
        "use_balanced": use_balanced,
        "normalization": "ICE-balanced" if use_balanced else "raw-count fallback",
    }
