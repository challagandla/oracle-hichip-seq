"""Run cooltools eigs-cis and normalise output to a stable TSV schema.

Different cooltools releases vary in whether -o is treated as a prefix or an
exact file. This wrapper always emits exactly snakemake.output.cis with columns
chrom, start, end, E1 so downstream scripts do not depend on cooltools naming.
"""
import subprocess
import sys
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from balance_utils import load_balance_report, resolution_balance  # noqa: E402
from utils import setup_logging, write_json  # noqa: E402


def _normalise_table(path: Path, out_path: Path) -> bool:
    try:
        df = pd.read_csv(path, sep="\t", comment="#")
    except Exception:
        return False
    lower = {c.lower(): c for c in df.columns}
    chrom = lower.get("chrom") or lower.get("chromosome")
    start = lower.get("start")
    end = lower.get("end")
    e1 = lower.get("e1") or lower.get("eig1") or lower.get("eigenvector")
    if not all([chrom, start, end, e1]):
        return False
    out = df[[chrom, start, end, e1]].rename(
        columns={chrom: "chrom", start: "start", end: "end", e1: "E1"}
    )
    out["chrom"] = out["chrom"].astype(str).str.strip()
    for column in ("start", "end", "E1"):
        out[column] = pd.to_numeric(out[column], errors="coerce")
    out = out.dropna(subset=["chrom", "start", "end", "E1"])
    out = out[
        (out["chrom"] != "")
        & (out["end"] > out["start"])
        & np.isfinite(out["E1"])
    ].copy()
    if out.empty:
        return False
    out[["start", "end"]] = out[["start", "end"]].astype(np.int64)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    return True


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    out_path = Path(snakemake.output.cis)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmpdir = out_path.parent / f".{snakemake.wildcards.sample}.eigs_tmp"
    shutil.rmtree(tmpdir, ignore_errors=True)
    res = int(snakemake.params.res)
    decision = resolution_balance(
        load_balance_report(snakemake.input.balance), res
    )
    status_payload = {
        "schema": "oracle-hichip-e1-status-v1",
        "sample": snakemake.wildcards.sample,
        "resolution_bp": res,
        "balance_status": decision["status"],
    }
    if not decision["use_balanced"]:
        # cooltools explicitly does not support eigs-cis on raw data. Preserve the
        # stable table schema without inventing zero-valued E1 measurements.
        pd.DataFrame(columns=["chrom", "start", "end", "E1"]).to_csv(
            out_path, sep="\t", index=False
        )
        status_payload.update({
            "status": "NOT_ASSESSED",
            "available": False,
            "normalization": None,
            "reason": "ICE balancing did not converge at the E1 resolution",
        })
        write_json(status_payload, snakemake.output.status)
        return

    tmpdir.mkdir(parents=True, exist_ok=True)
    prefix = tmpdir / "eigs"
    matrix = f"{snakemake.input.mcool}::resolutions/{res}"
    out_path.unlink(missing_ok=True)

    # No -p/--nproc: `cooltools expected-cis` and `cooltools insulation` take one,
    # but `eigs-cis` does not, and passing it aborts the command with
    # "Error: No such option '-p'" (cooltools 0.7.1).
    #
    # --view restricts the eigendecomposition to the assembled chromosomes. hg38
    # carries ~160 unplaced scaffolds; an A/B compartment call on a 60 kb contig is
    # meaningless even where it computes, and on a scaffold with no valid bins the
    # cooltools workers raise outright.
    # --phasing-track: the sign of an eigenvector is arbitrary, and cooltools solves
    # E1 one chromosome at a time, so without phasing "A" is positive on some
    # chromosomes and negative on others -- independently in each sample. Correlating
    # raw E1 across libraries can therefore average arbitrary sign flips toward zero
    # or produce negative replicate correlations. GC content is the standard phasing
    # track (the A compartment is generally GC-rich and gene-rich), and it fixes the
    # orientation consistently across samples.
    cmd = [
        "cooltools", "eigs-cis",
        "--view", str(snakemake.input.view),
        "--phasing-track", str(snakemake.input.gc),
        "--clr-weight-name", str(decision["weight_name"]),
        matrix, "-o", str(prefix),
    ]
    try:
        with open(snakemake.log[0], "a") as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=log)

        candidates = sorted(tmpdir.glob("*"))
        normalized = tmpdir / "normalized.tsv"
        for candidate in candidates:
            if candidate.is_file() and _normalise_table(candidate, normalized):
                normalized.replace(out_path)
                status_payload.update({
                    "status": "PASS",
                    "available": True,
                    "normalization": "ICE-balanced",
                })
                write_json(status_payload, snakemake.output.status)
                return
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    raise RuntimeError(f"cooltools eigs-cis completed but no usable E1 table was found for {snakemake.wildcards.sample}")


# Guarded so the module can be imported by the tests. Snakemake injects
# `snakemake` into the script's globals before executing it.
if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
