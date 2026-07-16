"""
Aggregate Peak Analysis (APA) on a HiChIP loop set.

For each loop, extract a square of the contact matrix centred on the loop and
average across loops. The centre pixel is then compared against two backgrounds:

  apa_vs_random_shift  centre / centre of random-shift controls where BOTH
                       anchors move by the same offset. Controls preserve
                       chromosome and distance, match sibling-anchor class and
                       caller-space marginal visibility, avoid blacklist loci,
                       and require a usable matrix window. This is a descriptive
                       effect size with a matched-loop bootstrap interval, not a
                       universal pass/fail statistic.
  apa_score            centre / two diagonal 3x3 corner neighbourhoods. Their
                       pixels span separation offsets within +/-2 bins of the
                       centre, so this is a near-distance corner ratio rather
                       than an exactly distance-matched control. Reported for
                       comparability with common APA summaries.

Contact frequency inside the window is a function of distance from the diagonal:
a pixel (i, j) lies at separation D + (j - i) * bin_size. So corners at different
(j - i) are not interchangeable, and averaging all four -- as this did -- puts a
corner that is 2 * win * bin_size CLOSER to the diagonal into the denominator.
"""
import logging
import sys
from pathlib import Path

import cooler
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from balance_utils import load_balance_report, resolution_balance  # noqa: E402
from utils import load_loops_bedpe, open_text_auto, setup_logging, write_json  # noqa: E402

log = logging.getLogger(__name__)

MAX_SHIFT_BP = 1_000_000


def _read_interval_index(
    path: str | Path,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Index BED intervals for O(log n) overlap checks."""
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    with open_text_auto(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 3:
                raise ValueError(f"invalid BED row in {path}: {line.rstrip()!r}")
            chrom, start, end = fields[:3]
            start_i, end_i = int(start), int(end)
            if start_i < 0 or end_i <= start_i:
                raise ValueError(f"invalid BED interval in {path}: {line.rstrip()!r}")
            by_chrom.setdefault(chrom, []).append((start_i, end_i))
    index: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for chrom, intervals in by_chrom.items():
        ordered = np.asarray(sorted(intervals), dtype=np.int64)
        index[chrom] = (ordered[:, 0], np.maximum.accumulate(ordered[:, 1]))
    return index


def _overlaps_interval(
    index: dict[str, tuple[np.ndarray, np.ndarray]],
    chrom: str,
    start: int,
    end: int,
) -> bool:
    if chrom not in index:
        return False
    starts, prefix_max_end = index[chrom]
    pos = int(np.searchsorted(starts, end, side="left") - 1)
    return pos >= 0 and int(prefix_max_end[pos]) > start


def _anchor_class(
    index: dict[str, tuple[np.ndarray, np.ndarray]],
    chrom: str,
    midpoint: int,
    bin_size: int,
) -> bool:
    start = midpoint // bin_size * bin_size
    return _overlaps_interval(index, chrom, start, start + bin_size)


def _caller_space_marginals(
    clr: cooler.Cooler,
    min_dist: int,
    max_dist: int,
    chunk_size: int = 1_000_000,
) -> np.ndarray:
    """Stream raw cis marginals over the FitHiChIP distance range."""
    bins = clr.bins()[:][["chrom", "start"]]
    n_bins = len(bins)
    starts = bins["start"].to_numpy(dtype=np.int64)
    chroms = bins["chrom"].astype(str).to_numpy()
    marginals = np.zeros(n_bins, dtype=float)
    selector = clr.pixels()
    nnz = int(clr.info["nnz"])
    if chunk_size < 1:
        raise ValueError("pixel chunk_size must be positive")
    for lo in range(0, nnz, int(chunk_size)):
        pixels = selector[lo:min(lo + int(chunk_size), nnz)][
            ["bin1_id", "bin2_id", "count"]
        ]
        bin1 = pixels["bin1_id"].to_numpy(dtype=np.int64)
        bin2 = pixels["bin2_id"].to_numpy(dtype=np.int64)
        distances = np.abs(starts[bin2] - starts[bin1])
        keep = (
            (chroms[bin1] == chroms[bin2])
            & (distances >= int(min_dist))
            & (distances <= int(max_dist))
        )
        values = pixels["count"].to_numpy(dtype=float)[keep]
        np.add.at(marginals, bin1[keep], values)
        np.add.at(marginals, bin2[keep], values)
    return marginals


def _bin_id(clr: cooler.Cooler, chrom: str, midpoint: int) -> int | None:
    try:
        lo, hi = clr.extent(chrom)
    except Exception:
        return None
    if midpoint < 0 or midpoint >= int(clr.chromsizes[chrom]):
        return None
    local = midpoint // int(clr.binsize)
    bin_id = int(lo) + int(local)
    return bin_id if int(lo) <= bin_id < int(hi) else None


def _marginally_matched(real: float, control: float, tolerance: float) -> bool:
    """Require nonzero caller-space visibility within a log2 fold tolerance."""
    if not np.isfinite(real) or not np.isfinite(control) or real <= 0 or control <= 0:
        return False
    return abs(float(np.log2(control / real))) <= float(tolerance)


def _matched_ratio_bootstrap(
    real_centres: np.ndarray,
    control_centres: np.ndarray,
    rng: np.random.Generator,
    n_boot: int = 500,
) -> tuple[float, float, float]:
    """Ratio and percentile CI from matched loop/control observations."""
    real = np.asarray(real_centres, dtype=float)
    control = np.asarray(control_centres, dtype=float)
    if real.shape != control.shape or real.ndim != 1 or real.size == 0:
        raise ValueError("APA bootstrap requires non-empty matched 1-D arrays")
    finite = np.isfinite(real) & np.isfinite(control)
    real = real[finite]
    control = control[finite]
    if real.size == 0:
        raise ValueError("APA bootstrap has no finite matched observations")
    point = float(real.mean() / max(control.mean(), 1e-9))
    boot = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        indices = rng.integers(0, real.size, size=real.size)
        boot[i] = real[indices].mean() / max(control[indices].mean(), 1e-9)
    low, high = np.quantile(boot, [0.025, 0.975])
    return point, float(low), float(high)


def _extract_square(clr: cooler.Cooler, chrom: str, mid1: int, mid2: int,
                    win_bp: int, balanced: bool = True) -> np.ndarray | None:
    """Extract a (2W+1)×(2W+1) window of contact frequency around (mid1,mid2)."""
    try:
        region1 = (chrom, max(0, mid1 - win_bp), mid1 + win_bp)
        region2 = (chrom, max(0, mid2 - win_bp), mid2 + win_bp)
        mat = clr.matrix(balance=balanced).fetch(region1, region2)
    except Exception:
        return None
    return mat


def _usable_square(square: np.ndarray | None, size: int, min_finite: float = 0.80) -> np.ndarray | None:
    """Reject truncated or mostly unmappable windows instead of zero-filling them."""
    if square is None:
        return None
    square = np.asarray(square, dtype=float)
    if square.shape != (size, size):
        return None
    finite = np.isfinite(square)
    centre = size // 2
    if finite.mean() < min_finite or not finite[centre, centre]:
        return None
    return square


def _not_assessed(snakemake, win: int, payload: dict, title: str) -> None:  # type: ignore[no-untyped-def]
    import matplotlib.pyplot as plt

    payload.update({"status": "NOT_ASSESSED", "pass": None})
    write_json(payload, snakemake.output.json)
    np.save(snakemake.output.npy, np.zeros((2 * win + 1, 2 * win + 1)))
    plt.figure()
    plt.title(title)
    plt.savefig(snakemake.output.png)
    plt.close()


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    import matplotlib.pyplot as plt

    setup_logging(snakemake.log[0])

    bin_sz = int(snakemake.params.bin_size)
    win = int(snakemake.params.window)            # in bins
    win_bp = win * bin_sz
    n_ctrl = int(snakemake.params.n_ctrl)
    balance = resolution_balance(
        load_balance_report(snakemake.input.balance), bin_sz
    )
    use_balanced = bool(balance["use_balanced"])
    analysis_context = {
        "normalization": balance["normalization"],
        "balance_status": balance["status"],
        "balance_converged": balance.get("converged"),
        "balance_resolution_bp": bin_sz,
        "candidate_source": (
            "loop calls from sibling donor(s), excluding the scored sample, "
            "restricted to sibling-only anchor support"
        ),
        "n_sibling_donor_callsets": int(snakemake.params.n_sibling_donor_callsets),
        "candidate_min_sample_support": int(snakemake.params.candidate_min_sample_support),
        "candidate_tolerance_bins": int(snakemake.params.candidate_tolerance_bins),
        "candidate_grid_bin_size_bp": int(
            snakemake.params.candidate_grid_bin_size_bp
        ),
        "candidate_reconciliation": (
            "q-filtered sibling reporting calls reconciled by reciprocal-anchor "
            "grid tolerance; separate from exact differential hypothesis pixels"
        ),
        "contact_map_held_out": bool(snakemake.params.contact_map_held_out),
        "anchors_exclude_scored_sample": bool(
            snakemake.params.anchors_exclude_scored_sample
        ),
        "primary_call_search_space": str(snakemake.params.primary_call_search_space),
        "candidate_loop_audit": str(snakemake.input.candidate_audit),
        "candidate_anchor_audit": str(snakemake.input.candidate_anchor_audit),
        "control_matching": {
            "same_chromosome_and_distance": True,
            "sibling_anchor_class": True,
            "blacklist_excluded": True,
            "usable_matrix_window": True,
            "caller_space_marginal_log2_tolerance": float(
                snakemake.params.control_marginal_log2_tolerance
            ),
            "visibility_distance_range_bp": [
                int(snakemake.params.visibility_min_dist),
                int(snakemake.params.visibility_max_dist),
            ],
        },
    }

    # The main diagonal must not enter the window, or it -- not the loop -- is what
    # the aggregate shows. A pixel at offset (dy, dx) from a loop of span D sits at
    # separation D + (dx - dy)*bin, and dx - dy ranges over +/- 2*win, so the diagonal
    # is inside the window for any loop with D <= 2*win*bin. A more permissive distance
    # floor can therefore turn the aggregate into a picture of the diagonal and inflate
    # apparent centre enrichment without measuring a loop-centred peak.
    #
    # Enforced here rather than trusted from config, so the geometry cannot be broken
    # again by editing a YAML value.
    min_dist_floor = (2 * win + 1) * bin_sz
    min_dist = int(snakemake.params.min_dist)
    if min_dist < min_dist_floor:
        log.warning(
            "apa.min_loop_dist=%d admits the main diagonal into a +/-%d-bin window; "
            "raising to %d (= (2*%d+1)*%d)",
            min_dist, win, min_dist_floor, win, bin_sz,
        )
        min_dist = min_dist_floor

    clr = cooler.Cooler(f"{snakemake.input.mcool}::resolutions/{bin_sz}")
    anchor_index = _read_interval_index(snakemake.input.candidate_anchors)
    blacklist_index = _read_interval_index(snakemake.input.blacklist)
    marginals = _caller_space_marginals(
        clr,
        int(snakemake.params.visibility_min_dist),
        int(snakemake.params.visibility_max_dist),
    )
    loops = load_loops_bedpe(snakemake.input.loops)

    n_called = len(loops)
    loops = loops[
        (loops.chrom1 == loops.chrom2) &
        ((loops.start2 - loops.start1).abs() >= min_dist)
    ].reset_index(drop=True)
    # Never silently. APA describes only the loops long enough to be measurable, and
    # how many were dropped to get there is part of the result.
    log.info(
        "APA on %d/%d loops (span >= %d bp; %d dropped as too short or trans)",
        len(loops), n_called, min_dist, n_called - len(loops),
    )

    # An aggregate over a handful of loops is not stable. One loop's centre divided by
    # a random-shift control that lands on empty matrix can produce an arbitrarily large
    # ratio. Below the floor the sample is NOT_ASSESSED: recorded, never passed.
    min_loops = int(snakemake.params.min_loops)
    if len(loops) < min_loops:
        log.warning(
            "only %d loops clear the distance floor (need %d); APA not assessed",
            len(loops), min_loops,
        )
        _not_assessed(snakemake, win, {
            **analysis_context,
            "sample": snakemake.wildcards.sample,
            "n_loops_called": int(n_called),
            "n_loops_used": int(len(loops)),
            "min_loops_for_apa": min_loops,
            "min_loop_dist_used": int(min_dist),
            "apa_score": None,
            "apa_vs_random_shift": None,
            "note": (
                f"Only {len(loops)} loops span >= {min_dist} bp; an aggregate over "
                f"fewer than {min_loops} loops reports noise, not enrichment."
            ),
        }, f"APA not assessed\n{len(loops)} held-out loops (need {min_loops})")
        return

    # Aggregate real loops. NaN bins are omitted per pixel, not converted to
    # contact zeros; zero-filling unmappable bins biases both centre and controls.
    real_records = []
    for loop_id, row in loops.iterrows():
        mid1 = (int(row.start1) + int(row.end1)) // 2
        mid2 = (int(row.start2) + int(row.end2)) // 2
        sq = _usable_square(
            _extract_square(
                clr, str(row.chrom1), mid1, mid2, win_bp,
                balanced=use_balanced,
            ),
            2 * win + 1,
        )
        if sq is not None:
            real_records.append((int(loop_id), row, sq))
    n_real_valid = len(real_records)
    if n_real_valid < min_loops:
        _not_assessed(snakemake, win, {
            **analysis_context,
            "sample": snakemake.wildcards.sample,
            "n_loops_called": int(n_called),
            "n_loops_candidates": int(len(loops)),
            "n_loops_used": int(n_real_valid),
            "min_loops_for_apa": min_loops,
            "min_loop_dist_used": int(min_dist),
            "apa_score": None,
            "apa_vs_random_shift": None,
            "note": "Too few held-out loop windows had adequate finite matrix coverage.",
        }, f"APA not assessed\n{n_real_valid} valid held-out windows")
        return

    # Generate controls only for loops with a valid observed window. Both anchors
    # move by the same offset, preserving chromosome and genomic distance. The
    # shifted pair must also preserve the two-anchor peak/non-peak class in the
    # sibling-only anchor universe, avoid blacklisted bins, and match each real
    # anchor's raw cis marginal in the FitHiChIP search range. This prevents a
    # peak-enriched HiChIP centre from being compared with an arbitrary invisible
    # non-peak locus. Repeated controls are reduced to one mean per loop so the
    # loop, rather than each random shift, remains the bootstrap unit.
    rng = np.random.default_rng(seed=42)
    control_centres: dict[int, list[float]] = {loop_id: [] for loop_id, _, _ in real_records}
    n_control_attempts = 0
    max_attempts = int(snakemake.params.max_control_attempts_per_draw)
    marginal_tolerance = float(snakemake.params.control_marginal_log2_tolerance)
    for loop_id, row, _ in real_records:
        chrom = str(row.chrom1)
        mid1 = (int(row.start1) + int(row.end1)) // 2
        mid2 = (int(row.start2) + int(row.end2)) // 2
        real_ids = (_bin_id(clr, chrom, mid1), _bin_id(clr, chrom, mid2))
        real_class = (
            _anchor_class(anchor_index, chrom, mid1, bin_sz),
            _anchor_class(anchor_index, chrom, mid2, bin_sz),
        )
        if None in real_ids or not any(real_class):
            continue
        real_visibility = (marginals[real_ids[0]], marginals[real_ids[1]])
        for _ in range(n_ctrl):
            for _attempt in range(max_attempts):
                n_control_attempts += 1
                # The offset must clear the window so the control does not overlap
                # the real APA centre.
                shift = int(rng.integers(win_bp + bin_sz, MAX_SHIFT_BP))
                if rng.random() < 0.5:
                    shift = -shift
                ctrl_mid1, ctrl_mid2 = mid1 + shift, mid2 + shift
                ctrl_ids = (
                    _bin_id(clr, chrom, ctrl_mid1),
                    _bin_id(clr, chrom, ctrl_mid2),
                )
                if None in ctrl_ids:
                    continue
                ctrl_class = (
                    _anchor_class(anchor_index, chrom, ctrl_mid1, bin_sz),
                    _anchor_class(anchor_index, chrom, ctrl_mid2, bin_sz),
                )
                if ctrl_class != real_class:
                    continue
                ctrl_bins = (
                    (ctrl_mid1 // bin_sz * bin_sz, ctrl_mid1 // bin_sz * bin_sz + bin_sz),
                    (ctrl_mid2 // bin_sz * bin_sz, ctrl_mid2 // bin_sz * bin_sz + bin_sz),
                )
                if any(
                    _overlaps_interval(blacklist_index, chrom, start, end)
                    for start, end in ctrl_bins
                ):
                    continue
                ctrl_visibility = (
                    marginals[ctrl_ids[0]], marginals[ctrl_ids[1]]
                )
                if not all(
                    _marginally_matched(real, control, marginal_tolerance)
                    for real, control in zip(real_visibility, ctrl_visibility)
                ):
                    continue
                sq = _usable_square(
                    _extract_square(
                        clr, chrom, ctrl_mid1, ctrl_mid2, win_bp,
                        balanced=use_balanced,
                    ),
                    2 * win + 1,
                )
                if sq is not None:
                    control_centres[loop_id].append(float(sq[win, win]))
                    break

    matched = [record for record in real_records if control_centres[record[0]]]
    n_used = len(matched)
    if n_used < min_loops:
        _not_assessed(snakemake, win, {
            **analysis_context,
            "sample": snakemake.wildcards.sample,
            "n_loops_called": int(n_called),
            "n_loops_with_valid_observed_windows": int(n_real_valid),
            "n_loops_used": int(n_used),
            "n_random_controls_requested": n_ctrl,
            "n_random_controls_usable": int(sum(map(len, control_centres.values()))),
            "n_random_control_attempts": int(n_control_attempts),
            "apa_score": None,
            "apa_vs_random_shift": None,
            "note": "Too few observed loops had at least one finite matched random-shift control.",
        }, f"APA not assessed\n{n_used} loops with matched controls")
        return

    real_squares = [record[2] for record in matched]
    with np.errstate(invalid="ignore"):
        agg = np.nanmean(np.stack(real_squares), axis=0)
    real_centres = np.asarray([record[2][win, win] for record in matched], dtype=float)
    per_loop_control = np.asarray(
        [np.mean(control_centres[record[0]]) for record in matched], dtype=float
    )
    apa_vs_ctrl, ci_low, ci_high = _matched_ratio_bootstrap(
        real_centres, per_loop_control, rng
    )

    centre = agg[win, win]

    # Background should remain near the centre pixel's genomic separation.
    # In this window a pixel (i, j) lies at separation D + (j - i) * bin_sz, so the
    # four corners span D ± 2 * win * bin_sz -- at win=20 and 10 kb bins that is a
    # 400 kb swing. The bottom-left corner is therefore ~400 kb CLOSER to the
    # diagonal than the loop, carries far more contacts for purely distance-decay
    # reasons, and dominates a four-corner mean: the old denominator was inflated
    # and the APA score correspondingly crushed. The two diagonal 3x3
    # neighbourhoods are centred on j-i=0; individual pixels retain only a
    # near-distance match within +/-2 bins. The random-shift control above is the
    # exact distance-matched effect used for primary reporting.
    background = np.concatenate([agg[:3, :3].ravel(), agg[-3:, -3:].ravel()])
    apa = float(centre / max(np.nanmean(background), 1e-9))

    # Plot
    fig, ax = plt.subplots(figsize=(4, 4))
    im = ax.imshow(np.log2(agg + 1), cmap="Reds", origin="lower",
                   extent=[-win, win, -win, win])
    ax.set_title(
        f"Contact-map-held-out APA {snakemake.wildcards.sample}\n"
        f"matched effect={apa_vs_ctrl:.2f} (95% bootstrap interval "
        f"{ci_low:.2f}-{ci_high:.2f}), n={n_used}\n"
        f"{balance['normalization']} (balance {balance['status']})"
    )
    ax.set_xlabel(f"bins ({bin_sz//1000} kb)")
    ax.set_ylabel("bins")
    fig.colorbar(im, ax=ax, label="log2(1+contacts)")
    fig.tight_layout()
    fig.savefig(snakemake.output.png, dpi=150)

    # The matched random-shift ratio is reported as a descriptive effect size, not
    # a universal pass/fail gate. Its interval quantifies loop-resampling uncertainty
    # within this dataset; it does not calibrate a cross-study quality threshold.
    np.save(snakemake.output.npy, agg)

    write_json({
        **analysis_context,
        "sample": snakemake.wildcards.sample,
        "n_loops_called": int(n_called),
        "n_loops_used": int(n_used),
        "n_random_controls_requested": n_ctrl,
        "n_random_controls_usable": int(sum(map(len, control_centres.values()))),
        "n_random_control_attempts": int(n_control_attempts),
        "median_controls_per_loop": float(
            np.median([len(control_centres[record[0]]) for record in matched])
        ),
        # The distance floor and the loss it causes travel with the score. An APA of
        # 3.5 over 14,000 loops and an APA of 3.5 over 60 are not the same statement,
        # and neither is comparable to one whose window straddled the diagonal.
        "min_loop_dist_used": int(min_dist),
        "apa_window_bins": int(win),
        "apa_score": apa,
        "apa_corner_background": (
            "two diagonal 3x3 corner neighbourhoods; separation offset within "
            "+/-2 bins (near-distance, not exact distance matching)"
        ),
        "apa_vs_random_shift": apa_vs_ctrl,
        "apa_vs_random_shift_ci95": [ci_low, ci_high],
        "status": "DESCRIPTIVE",
        "pass": None,
        "note": (
            "APA is a descriptive matched effect size and is never a hard sample "
            "pass/fail gate. ICE balancing did not converge, so this value uses an "
            "explicitly labelled raw-count fallback."
            if not use_balanced else
            "APA is a descriptive matched effect size and is never a hard sample "
            "pass/fail gate."
        ),
    }, snakemake.output.json)


# Guarded so the module can be imported by the tests. Snakemake injects
# `snakemake` into the script's globals before executing it.
if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
