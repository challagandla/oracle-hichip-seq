"""
Stratum-adjusted Pearson correlation (HiCRep) between biological replicates.

Outputs a four-state QC result: PASS, FAIL, DISCORDANT, or NOT_ASSESSED.
Single-replicate samples are not treated as true passes. For groups with more
than two replicates, a favourable best match cannot hide a weak pair: every
depth-qualified pair must clear the threshold for PASS.

hicrepSCC returns a numpy array indexed by chromosome, pre-filled with the
sentinel -2.0 (`scc = np.full(len(chrNames), -2.0)`), not a dict and not NaN. So
chromosomes it could not score keep -2.0, and averaging without masking them
drags SCC toward -2. Chromosomes are excluded via the function's own excludeChr
argument rather than by matching substrings of the returned keys.
"""
import logging
import hashlib
import sys
from itertools import combinations
from pathlib import Path

import cooler
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from utils import setup_logging, write_json  # noqa: E402

SCC_SENTINEL = -2.0

log = logging.getLogger(__name__)


def _scc_chroms(view_path: str, clr_chroms: set[str]) -> list[str]:
    """Autosomes present in the cooler, taken from the shared main-chromosome view.

    hicrepSCC does not skip a chromosome it cannot score -- it asserts:

        AssertionError: Contact matrix 1 of chromosome GL000208.1 is empty

    and hg38 has ~160 unplaced scaffolds that are empty at 25 kb, so passing the
    cooler's full chromosome list kills the rule. Naming the chromosomes explicitly
    is also the correct thing statistically: SCC is a stratum-adjusted correlation
    over a distance-decay profile, which a 60 kb contig does not have.

    chrX is dropped even though it is in the view. An autosomal SCC contract avoids
    sex-chromosome ploidy and mappability effects, remains comparable when users
    replace the bundled male cohort, and matches the depth denominator exactly.
    """
    keep = []
    for line in Path(view_path).read_text().splitlines():
        if not line.strip():
            continue
        c = line.split("\t")[0]
        if c in clr_chroms and c != "chrX":
            keep.append(c)
    return keep


def _selected_cis_contacts(
    clr: cooler.Cooler, chromosomes: list[str], max_dist: int
) -> int:
    """Contacts in the exact diagonal/chromosome population scored by HiCRep.

    hicrepSCC removes the main diagonal and diagonals beyond ``dBPMax`` before
    stochastic depth matching. Counting all cis pixels would let long-range or
    diagonal contacts clear a depth floor even though they never enter SCC.
    """
    pixels = clr.matrix(balance=False, as_pixels=True, join=False)
    max_bin_offset = int(max_dist) // int(clr.binsize)
    total = 0
    for chrom in chromosomes:
        table = pixels.fetch(chrom)
        offsets = table["bin2_id"].to_numpy() - table["bin1_id"].to_numpy()
        keep = (offsets > 0) & (offsets <= max_bin_offset)
        total += int(table.loc[keep, "count"].sum())
    return total


def _cis_contacts(
    path: str, bin_sz: int, chromosomes: list[str], max_dist: int
) -> int:
    clr = cooler.Cooler(f"{path}::resolutions/{bin_sz}")
    return _selected_cis_contacts(clr, chromosomes, max_dist)


def _pair_seed(a: str, b: str, bin_sz: int, max_dist: int, h: int) -> int:
    """Stable NumPy seed for HiCRep's stochastic depth matching."""
    identity = "|".join([*sorted((Path(a).stem, Path(b).stem)), str(bin_sz), str(max_dist), str(h)])
    return int(hashlib.sha256(identity.encode()).hexdigest()[:8], 16)


def _classify_sccs(
    values: list[float], threshold: float
) -> tuple[str, bool | None]:
    """Classify a set of depth-qualified replicate comparisons.

    ``best_scc`` is intentionally not used here. Taking the maximum is a
    replicate-count-dependent selection statistic: one strong sibling match can
    conceal another discordant pair. Mixed evidence is therefore explicit rather
    than forced into either PASS or FAIL.
    """
    if not values:
        return "NOT_ASSESSED", None
    passes = [value >= threshold for value in values]
    if all(passes):
        return "PASS", True
    if not any(passes):
        return "FAIL", False
    return "DISCORDANT", None


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    from hicrep import hicrepSCC
    from hicrep.utils import readMcool
    mcools = list(snakemake.input.mcools)
    bin_sz = int(snakemake.params.bin)
    max_dist = int(snakemake.params.maxd)
    h = int(snakemake.params.h)
    threshold = float(snakemake.params.threshold)
    min_contacts = int(snakemake.params.min_contacts)

    clr_chroms = set(cooler.Cooler(f"{mcools[0]}::resolutions/{bin_sz}").chromnames)
    scc_chroms = _scc_chroms(snakemake.input.view, clr_chroms)
    log.info("scoring SCC over %d chromosomes", len(scc_chroms))
    depth = {
        p: _cis_contacts(p, bin_sz, scc_chroms, max_dist) for p in mcools
    }

    this_sample = snakemake.wildcards.sample

    result = {
        "sample": this_sample,
        "bin_size": bin_sz,
        "max_dist": max_dist,
        "h": h,
        "n_replicates": len(mcools),
        "contacts": {Path(p).stem: depth[p] for p in mcools},
        "contact_depth_population": (
            "stored off-diagonal cis contacts on HiCRep-scored autosomes at "
            f"1..{max_dist // bin_sz} bins (<= {max_dist} bp), exactly matching "
            "the population retained before hicrepSCC depth matching"
        ),
        "min_contacts_for_scc": min_contacts,
        "pairwise_scc": [],
        "mean_scc": None,
        "min_scc": None,
        "best_scc": None,
        "n_qualified_pairs": 0,
        "group_median_scc": None,
        "group_n_qualified_pairs": 0,
        "group_status": "NOT_ASSESSED",
        "group_pass": None,
        "threshold": threshold,
        "status": "NOT_ASSESSED",
        "pass": None,
    }

    if len(mcools) < 2:
        result["note"] = "Only one replicate; HiCRep was not assessed."
        write_json(result, snakemake.output.json)
        return

    sccs = []
    group_sccs = []
    for a, b in combinations(mcools, 2):
        # SCC is depth-sensitive because bDownSample reduces the deeper library to
        # the shallower contact count. Below the configured floor, report the value
        # as depth-confounded rather than using it for PASS/FAIL.
        shallow = min(depth[a], depth[b])
        confounded = shallow < min_contacts

        cool_a, _ = readMcool(a, bin_sz)
        cool_b, _ = readMcool(b, bin_sz)
        # bDownSample=True: HiCRep's SCC is sensitive to sequencing depth, and
        # the deeper matrix is downsampled to the shallower contact count.
        seed = _pair_seed(a, b, bin_sz, max_dist, h)
        np.random.seed(seed)
        scc = np.asarray(
            hicrepSCC(cool_a, cool_b, h, max_dist, True, chrNames=scc_chroms),
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
            "status": (
                "DEPTH_CONFOUNDED"
                if confounded
                else ("PASS" if mean_scc >= threshold else "FAIL")
            ),
            "downsample_seed": seed,
        })
        if confounded:
            log.warning(
                "%s vs %s: SCC=%.3f but the shallower library holds %d contacts "
                "(floor %d) -- reporting as depth-confounded, not as concordance",
                Path(a).stem, Path(b).stem, mean_scc, shallow, min_contacts,
            )
            continue
        group_sccs.append(mean_scc)
        # A per-sample report uses only pairs that contain that sample; group-level
        # pairs are retained separately for cohort context.
        if this_sample in (Path(a).stem, Path(b).stem):
            sccs.append(mean_scc)

    result["group_n_qualified_pairs"] = len(group_sccs)
    result["group_status"], result["group_pass"] = _classify_sccs(
        group_sccs, threshold
    )
    if group_sccs:
        result["group_median_scc"] = float(np.median(group_sccs))

    if not sccs:
        result["note"] = (
            "No replicate pair involving this library both cleared the depth floor "
            f"of {min_contacts} contacts at {bin_sz} bp and produced a usable SCC; "
            "replicate concordance was not assessed."
        )
        write_json(result, snakemake.output.json)
        return

    result["n_qualified_pairs"] = len(sccs)
    result["mean_scc"] = float(np.mean(sccs))
    result["min_scc"] = float(np.min(sccs))
    result["best_scc"] = float(np.max(sccs))
    result["status"], result["pass"] = _classify_sccs(sccs, threshold)
    if result["status"] == "DISCORDANT":
        result["note"] = (
            "Depth-qualified replicate pairs disagree about the SCC threshold; "
            "inspect every pair rather than selecting the best match."
        )
    write_json(result, snakemake.output.json)


# Guarded so the module can be imported by the tests. Snakemake injects
# `snakemake` into the script's globals before executing it.
if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
