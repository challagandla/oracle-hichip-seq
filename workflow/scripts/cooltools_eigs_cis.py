"""Run cooltools eigs-cis and normalise output to a stable TSV schema.

Different cooltools releases vary in whether -o is treated as a prefix or an
exact file. This wrapper always emits exactly snakemake.output.cis with columns
chrom, start, end, E1 so downstream scripts do not depend on cooltools naming.
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from utils import setup_logging  # noqa: E402


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
    out = df[[chrom, start, end, e1]].rename(columns={chrom: "chrom", start: "start", end: "end", e1: "E1"})
    out = out.dropna(subset=["chrom", "start", "end", "E1"])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    return True


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    out_path = Path(snakemake.output.cis)
    tmpdir = out_path.parent / f".{snakemake.wildcards.sample}.eigs_tmp"
    tmpdir.mkdir(parents=True, exist_ok=True)
    prefix = tmpdir / "eigs"
    res = int(snakemake.params.res)
    matrix = f"{snakemake.input.mcool}::resolutions/{res}"

    # No -p/--nproc: `cooltools expected-cis` and `cooltools insulation` take one,
    # but `eigs-cis` does not, and passing it aborts the command with
    # "Error: No such option '-p'" (cooltools 0.7.1).
    #
    # --view restricts the eigendecomposition to the assembled chromosomes. hg38
    # carries ~160 unplaced scaffolds; an A/B compartment call on a 60 kb contig is
    # meaningless even where it computes, and on a scaffold with no valid bins the
    # cooltools workers raise outright.
    cmd = [
        "cooltools", "eigs-cis",
        "--view", str(snakemake.input.view),
        matrix, "-o", str(prefix),
    ]
    with open(snakemake.log[0], "a") as log:
        subprocess.run(cmd, check=True, stdout=log, stderr=log)

    candidates = []
    if out_path.exists():
        candidates.append(out_path)
    candidates.extend(sorted(tmpdir.glob("*")))
    candidates.extend(sorted(out_path.parent.glob(f"{prefix.name}*")))

    for candidate in candidates:
        if candidate.is_file() and _normalise_table(candidate, out_path):
            return

    raise RuntimeError(f"cooltools eigs-cis completed but no usable E1 table was found for {snakemake.wildcards.sample}")


# Guarded so the module can be imported by the tests. Snakemake injects
# `snakemake` into the script's globals before executing it.
if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
