"""
Stratum-adjusted Pearson correlation (HiCRep) between biological replicates.

Outputs a three-state QC result: PASS, FAIL, or NOT_ASSESSED. Single-replicate
samples are not treated as true passes.

hicrepSCC returns a numpy array indexed by chromosome, pre-filled with the
sentinel -2.0 (`scc = np.full(len(chrNames), -2.0)`), not a dict and not NaN. So
chromosomes it could not score keep -2.0, and averaging without masking them
drags SCC toward -2. Chromosomes are excluded via the function's own excludeChr
argument rather than by matching substrings of the returned keys.
"""
import logging
import sys
from itertools import combinations
from pathlib import Path

import cooler
import numpy as np
from hicrep import hicrepSCC
from hicrep.utils import readMcool

sys.path.insert(0, str(Path(__file__).parent))
from utils import setup_logging, write_json  # noqa: E402

# Sex chromosomes differ in copy number between donors and chrM has no meaningful
# contact structure; either would corrupt a between-donor correlation.
EXCLUDE_CHROMS = {"chrX", "chrY", "chrM", "X", "Y", "M", "MT"}
SCC_SENTINEL = -2.0

log = logging.getLogger(__name__)


def _cis_contacts(path: str, bin_sz: int) -> int:
    """Total contacts in the matrix HiCRep will actually read."""
    return int(cooler.Cooler(f"{path}::resolutions/{bin_sz}").info["sum"])


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    mcools = list(snakemake.input.mcools)
    bin_sz = int(snakemake.params.bin)
    max_dist = int(snakemake.params.maxd)
    h = int(snakemake.params.h)
    threshold = float(snakemake.config["hicrep"]["threshold_pass"])
    min_contacts = int(snakemake.config["hicrep"]["min_contacts_for_scc"])

    depth = {p: _cis_contacts(p, bin_sz) for p in mcools}

    result = {
        "sample": snakemake.wildcards.sample,
        "bin_size": bin_sz,
        "max_dist": max_dist,
        "h": h,
        "n_replicates": len(mcools),
        "contacts": {Path(p).stem: depth[p] for p in mcools},
        "min_contacts_for_scc": min_contacts,
        "pairwise_scc": [],
        "mean_scc": None,
        "threshold": threshold,
        "status": "NOT_ASSESSED",
        "pass": None,
    }

    if len(mcools) < 2:
        result["note"] = "Only one replicate; HiCRep was not assessed."
        write_json(result, snakemake.output.json)
        return

    sccs = []
    for a, b in combinations(mcools, 2):
        # SCC is dominated by the shallower library, because bDownSample crushes the
        # deeper matrix to the shallower one's contact count. Measured on this
        # cohort at 25 kb: two libraries from DIFFERENT cell types (90.8M and 47.6M
        # contacts) score 0.743, while two genuine replicates score 0.221 when one
        # of them holds 3.0M. Below the floor the number reports depth, not
        # concordance, so it is recorded but never allowed to decide PASS/FAIL.
        shallow = min(depth[a], depth[b])
        confounded = shallow < min_contacts

        cool_a, _ = readMcool(a, bin_sz)
        cool_b, _ = readMcool(b, bin_sz)
        # bDownSample=True: HiCRep's SCC is sensitive to sequencing depth, and
        # these libraries differ by an order of magnitude, so the deeper matrix is
        # downsampled to the shallower one's contact count before comparison.
        # It was previously being passed `bin_sz` positionally into this slot.
        scc = np.asarray(
            hicrepSCC(cool_a, cool_b, h, max_dist, True, excludeChr=EXCLUDE_CHROMS),
            dtype=float,
        )
        scored = scc[(scc > SCC_SENTINEL) & np.isfinite(scc)]
        if scored.size == 0:
            log.warning("HiCRep scored no chromosome for %s vs %s", a, b)
            continue
        mean_scc = float(scored.mean())
        result["pairwise_scc"].append({
            "a": Path(a).stem, "b": Path(b).stem,
            "scc": mean_scc, "n_chroms_scored": int(scored.size),
            "min_contacts": int(shallow), "depth_confounded": bool(confounded),
        })
        if confounded:
            log.warning(
                "%s vs %s: SCC=%.3f but the shallower library holds %d contacts "
                "(floor %d) -- reporting as depth-confounded, not as concordance",
                Path(a).stem, Path(b).stem, mean_scc, shallow, min_contacts,
            )
            continue
        sccs.append(mean_scc)

    if not sccs:
        result["note"] = (
            "No replicate pair cleared the depth floor of "
            f"{min_contacts} contacts at {bin_sz} bp; SCC would report sequencing "
            "depth rather than replicate concordance, so it was not assessed."
        )
        write_json(result, snakemake.output.json)
        return

    result["mean_scc"] = float(np.mean(sccs))
    result["pass"] = result["mean_scc"] >= threshold
    result["status"] = "PASS" if result["pass"] else "FAIL"
    write_json(result, snakemake.output.json)


# Guarded so the module can be imported by the tests. Snakemake injects
# `snakemake` into the script's globals before executing it.
if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
