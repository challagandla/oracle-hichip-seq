"""
Cohort-level publication figures.

The per-sample PNGs the other stages emit are diagnostics. These are the figures
that carry the argument, so they are assembled across the whole cohort, written as
vector PDF alongside 400 dpi PNG. Scientifically valid empty results render an
explicit "no data" panel; missing or malformed required inputs fail the rule.
That distinction prevents a truncated report from being mistaken for a negative
biological result.

Nothing here silently hides a bad library. Shallow libraries stay on the plots and
are marked, because the reader's first question about a HiChIP cohort is whether
the contrast is driven by biology or by depth.
"""
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
except ModuleNotFoundError:  # Pure reporting-contract helpers remain importable.
    matplotlib = None
    plt = None
    Line2D = None

sys.path.insert(0, str(Path(__file__).parent))
from utils import setup_logging  # noqa: E402

log = logging.getLogger(__name__)

# ---------------------------------------------------------------- style
PLOT_STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 8,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
    "savefig.bbox": "tight",
}
if plt is not None:
    plt.rcParams.update(PLOT_STYLE)

PALETTE = {"Naive": "#4C72B0", "Th17": "#C44E52", "Treg": "#55A868"}
GREY = "#9e9e9e"


def _require_plotting() -> None:
    if plt is None or Line2D is None:
        raise RuntimeError(
            "matplotlib is required to render publication figures; use the "
            "workflow figures environment"
        )


def _colour(cell_type: str) -> str:
    return PALETTE.get(cell_type, GREY)


def _panel_label(ax, letter: str) -> None:
    ax.text(-0.22, 1.06, letter, transform=ax.transAxes,
            fontsize=10, fontweight="bold", va="bottom", ha="left")


def _short(sample: str) -> str:
    """Naive_H3K27ac_rep1 -> 'Naive K27ac r1'.

    Full sample ids are long enough that, rotated 90 degrees, they take more of the
    panel than the data does.
    """
    s = str(sample).replace("H3K27ac", "K27ac").replace("H3K4me3", "K4me3")
    s = s.replace("H3K4me1", "K4me1").replace("H3K27me3", "K27me3")
    return s.replace("_rep", " r").replace("_", " ")


def _label_axis(ax, samples, axis: str = "x") -> None:
    labels = [_short(s) for s in samples]
    if axis == "x":
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=90, fontsize=5)
    else:
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=5)


def primary_reporting_libraries(lib: pd.DataFrame) -> pd.DataFrame:
    """Return the libraries eligible for headline biological reporting."""
    if "report_role" not in lib.columns:
        raise ValueError("library table lacks report_role")
    return lib[lib["report_role"] == "primary"].copy()


def qc_gate_values(thresholds: dict) -> tuple[float, float, float]:
    """Return configured yield, duplicate, and cis-percent gates."""
    required = {"valid_pair_yield_pct", "duplicate_pct_max", "cis_fraction_min"}
    missing = sorted(required - set(thresholds))
    if missing:
        raise ValueError(f"configured QC thresholds are missing: {missing}")
    return (
        float(thresholds["valid_pair_yield_pct"]),
        float(thresholds["duplicate_pct_max"]),
        100.0 * float(thresholds["cis_fraction_min"]),
    )


def qc_gate_failures(
    values: pd.Series, *, minimum: float | None = None,
    maximum: float | None = None,
) -> pd.Series:
    """Return which assessed values fail one directional configured gate."""
    if (minimum is None) == (maximum is None):
        raise ValueError("provide exactly one of minimum or maximum")
    numeric = pd.to_numeric(values, errors="coerce")
    if minimum is not None:
        return numeric.notna() & (numeric < float(minimum))
    return numeric.notna() & (numeric > float(maximum))


def _empty(ax, msg: str) -> None:
    ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=6,
            color=GREY, transform=ax.transAxes, wrap=True)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def _save(fig, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=400)
    plt.close(fig)
    log.info("wrote %s.{pdf,png}", stem)


# ---------------------------------------------------------------- parsing
def _required_file(path: str | Path, label: str, *, allow_empty: bool = False) -> Path:
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"required {label} is missing: {source}")
    if not allow_empty and source.stat().st_size == 0:
        raise ValueError(f"required {label} is empty: {source}")
    return source


def _required_json(
    path: str | Path, label: str, required_keys: set[str] | tuple[str, ...]
) -> dict:
    source = _required_file(path, label)
    try:
        value = json.loads(source.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"required {label} is not valid JSON: {source}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"required {label} must be a JSON object: {source}")
    missing = sorted(set(required_keys) - set(value))
    if missing:
        raise ValueError(f"required {label} lacks keys {missing}: {source}")
    return value


def _required_tsv(
    path: str | Path,
    label: str,
    required_columns: set[str] | tuple[str, ...],
    *,
    allow_empty_rows: bool = True,
) -> pd.DataFrame:
    source = _required_file(path, label)
    try:
        frame = pd.read_csv(source, sep="\t")
    except Exception as exc:
        raise ValueError(f"required {label} is not a valid TSV: {source}") from exc
    missing = sorted(set(required_columns) - set(frame.columns))
    if missing:
        raise ValueError(f"required {label} lacks columns {missing}: {source}")
    if frame.empty and not allow_empty_rows:
        raise ValueError(f"required {label} has no data rows: {source}")
    return frame


def _validate_loop_coordinates(frame: pd.DataFrame, label: str) -> pd.DataFrame:
    """Validate canonical BEDPE coordinates without coercing or dropping rows."""
    required = {"chrom1", "start1", "end1", "chrom2", "start2", "end2"}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"required {label} lacks columns {missing}")
    frame = frame.copy()
    coordinates = ("start1", "end1", "start2", "end2")
    for column in coordinates:
        try:
            values = pd.to_numeric(frame[column], errors="raise")
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"required {label} has nonnumeric {column} coordinates"
            ) from exc
        if values.isna().any() or not np.equal(values, np.floor(values)).all():
            raise ValueError(f"required {label} has non-integer {column} coordinates")
        frame[column] = values.astype(np.int64)
    if (
        (frame[["start1", "start2"]] < 0).any().any()
        or (frame["end1"] <= frame["start1"]).any()
        or (frame["end2"] <= frame["start2"]).any()
    ):
        raise ValueError(f"required {label} has invalid BEDPE intervals")
    for column in ("chrom1", "chrom2"):
        if frame[column].isna().any() or (frame[column].astype(str).str.strip() == "").any():
            raise ValueError(f"required {label} has empty {column} values")
        frame[column] = frame[column].astype(str)
    return frame


def _required_pipeline_loops(path: str | Path, label: str) -> pd.DataFrame:
    """Read canonical pipeline BEDPE without dropping malformed coordinates."""
    frame = _required_tsv(
        path,
        label,
        {"chrom1", "start1", "end1", "chrom2", "start2", "end2"},
    )
    return _validate_loop_coordinates(frame, label)


def validate_apa_matrix_contract(metadata: dict, matrix: np.ndarray, sample: str) -> str:
    """Cross-check the APA categorical state against its numeric matrix."""
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or matrix.shape[0] % 2 != 1:
        raise ValueError(
            f"APA matrix for {sample} must be an odd square 2-D array; got {matrix.shape}"
        )
    status = str(metadata.get("status"))
    if status == "NOT_ASSESSED":
        if matrix.any():
            raise ValueError(
                f"APA metadata says NOT_ASSESSED but matrix is nonzero for {sample}"
            )
        return status
    if status != "DESCRIPTIVE":
        raise ValueError(f"unknown APA status {status!r} for {sample}")
    if not matrix.any():
        raise ValueError(
            f"APA metadata says DESCRIPTIVE but matrix is all zero for {sample}"
        )
    return status


def _pairtools_stats(path: Path, required: set[str] | None = None) -> dict[str, float]:
    out: dict[str, float] = {}
    source = _required_file(path, "pairtools statistics")
    for line in source.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            try:
                out[parts[0]] = float(parts[1])
            except ValueError:
                pass
    missing = sorted((required or set()) - set(out))
    if missing:
        raise ValueError(f"pairtools statistics lack numeric keys {missing}: {source}")
    return out


def _library_table(
    samples: pd.DataFrame,
    results: Path,
    min_contacts: int,
    demonstration_samples: set[str] | None = None,
) -> pd.DataFrame:
    demonstration_samples = demonstration_samples or set()
    rows = []
    for sid, meta in samples.iterrows():
        dedup = _pairtools_stats(
            results / f"qc/pairtools/{sid}.dedup.stats.txt",
            {"total", "total_mapped", "total_dups"},
        )
        selected = _pairtools_stats(
            results / f"qc/pairtools/{sid}.pairs.stats.txt", {"total", "cis"}
        )
        pairtools_input_pairs = dedup.get("total", np.nan)
        unique_uu = selected.get("total", selected.get("total_nodups", np.nan))
        dup_fraction = dedup.get("summary/frac_dups", np.nan)
        if np.isnan(dup_fraction):
            mapped = dedup.get("total_mapped", np.nan)
            dups = dedup.get("total_dups", np.nan)
            dup_fraction = dups / mapped if mapped and not np.isnan(mapped) else np.nan
        cis_fraction = selected.get("summary/frac_cis", np.nan)
        if np.isnan(cis_fraction):
            cis = selected.get("cis", np.nan)
            cis_fraction = cis / unique_uu if unique_uu and not np.isnan(unique_uu) else np.nan

        loops_f = results / f"loops/{sid}/{sid}.interactions_FitHiC_{Q_LABEL}.bed"
        n_loops = len(
            _required_pipeline_loops(loops_f, f"FitHiChIP loops for {sid}")
        )

        apa_f = results / f"qc/apa/{sid}.apa.json"
        apa = _required_json(
            apa_f,
            f"APA metadata for {sid}",
            {"status", "normalization", "apa_vs_random_shift"},
        )
        balance_f = results / f"qc/balance/{sid}.balance.json"
        balance = _required_json(balance_f, f"balance QC for {sid}", {"status"})
        loop_qc_f = results / f"qc/loop_qc/{sid}.json"
        loop_qc = _required_json(
            loop_qc_f,
            f"loop QC for {sid}",
            {
                "raw_input_pairs", "post_trim_pairs", "valid_pair_yield_pct",
                "post_trim_valid_pair_yield_pct", "overall_status",
            },
        )
        raw_input_pairs = loop_qc.get("raw_input_pairs", np.nan)
        post_trim_pairs = loop_qc.get("post_trim_pairs", np.nan)
        contact_depth_f = results / f"qc/contact_depth/{sid}.json"
        contact_depth = _required_json(
            contact_depth_f,
            f"contact-depth QC for {sid}",
            {
                "primary_cis_offdiagonal_contacts",
                "fithichip_distance_range_contacts",
            },
        )
        hicrep_f = results / f"qc/hicrep/{sid}.hicrep.json"
        hicrep = _required_json(
            hicrep_f, f"HiCRep QC for {sid}", {"status", "contacts", "pairwise_scc"}
        )
        hicrep_contacts = (hicrep.get("contacts") or {}).get(sid)

        rows.append({
            "sample": sid,
            "report_role": (
                "demonstration" if sid in demonstration_samples else "primary"
            ),
            "cell_type": meta.get("cell_type", "?"),
            "mark": meta.get("mark", "?"),
            "raw_input_pairs": raw_input_pairs,
            "post_trim_pairs": post_trim_pairs,
            "pairtools_input_pairs": pairtools_input_pairs,
            "unique_pairs": unique_uu,
            "valid_yield_pct": loop_qc.get("valid_pair_yield_pct", np.nan),
            "post_trim_valid_yield_pct": loop_qc.get(
                "post_trim_valid_pair_yield_pct", np.nan
            ),
            "dup_pct": 100.0 * dup_fraction,
            "cis_pct": 100.0 * cis_fraction,
            "n_loops": n_loops,
            "n_consensus_peaks": loop_qc.get("n_consensus_peaks"),
            "loop_search_contacts": contact_depth.get(
                "fithichip_distance_range_contacts"
            ),
            "stripe_search_contacts": contact_depth.get(
                "primary_cis_offdiagonal_contacts"
            ),
            "apa": apa.get("apa_vs_random_shift"),
            "apa_status": apa.get("status", "NOT_ASSESSED"),
            "apa_normalization": apa.get("normalization"),
            "balance_status": balance.get("status", "NOT_ASSESSED"),
            "overall_qc_status": loop_qc.get("overall_status", "NOT_ASSESSED"),
            "hicrep_contacts": hicrep_contacts,
            "hicrep_depth_adequate": bool(
                hicrep_contacts is not None and hicrep_contacts >= min_contacts
            ),
        })
    return pd.DataFrame(rows).set_index("sample")


# ---------------------------------------------------------------- figure 1
def figure1_library_qc(
    lib: pd.DataFrame,
    results: Path,
    out: Path,
    min_contacts: int,
    thresholds: dict,
) -> None:
    yield_gate, duplicate_gate, cis_gate = qc_gate_values(thresholds)

    fig = plt.figure(figsize=(7.6, 6.2))
    gs = fig.add_gridspec(2, 3, hspace=1.05, wspace=0.62)

    # (a) Raw valid-pair yield: final UU contacts / all sequenced read pairs.
    ax = fig.add_subplot(gs[0, 0])
    _panel_label(ax, "a")
    d = lib.dropna(subset=["valid_yield_pct"]).sort_values("valid_yield_pct")
    if len(d):
        bars = ax.barh(
            range(len(d)), d["valid_yield_pct"],
            color=[_colour(c) for c in d["cell_type"]],
            edgecolor="black", linewidth=0.3,
        )
        failed = qc_gate_failures(d["valid_yield_pct"], minimum=yield_gate)
        for bar, is_failed in zip(bars, failed):
            if is_failed:
                bar.set_hatch("///")
                bar.set_edgecolor("#b00020")
        ax.axvline(yield_gate, ls="--", lw=0.7, c="#b00020")
        _label_axis(ax, d.index, axis="y")
        ax.set_xlim(0, max(100, float(d["valid_yield_pct"].max()) * 1.05))
        ax.set_xlabel("raw valid-pair yield (%)")
        ax.set_title("End-to-end contact yield", loc="left")
    else:
        _empty(ax, "no valid-pair yield records")

    # (b) HiCRep-scored cis contacts, with the matching depth floor drawn on.
    ax = fig.add_subplot(gs[0, 1])
    _panel_label(ax, "b")
    d = lib.dropna(subset=["hicrep_contacts"]).sort_values(
        "hicrep_contacts", ascending=False
    )
    if len(d):
        colours = [_colour(c) for c in d["cell_type"]]
        bars = ax.bar(range(len(d)), d["hicrep_contacts"] / 1e6, color=colours,
                      edgecolor="black", linewidth=0.3)
        for i, (_, r) in enumerate(d.iterrows()):
            if not r["hicrep_depth_adequate"]:
                bars[i].set_hatch("///")
                bars[i].set_edgecolor("#b00020")
        ax.axhline(min_contacts / 1e6, ls="--", lw=0.7, c="#b00020")
        # Annotate above the axes, not inside them: the line sits low on a log axis
        # and the label lands on top of the tick labels.
        ax.text(0.99, 1.01, f"{min_contacts/1e6:.0f}M HiCRep depth floor",
                transform=ax.transAxes, fontsize=5, c="#b00020", ha="right", va="bottom")
        _label_axis(ax, d.index, axis="x")
        ax.set_ylabel("HiCRep-scored cis contacts (millions)")
        ax.set_yscale("log")
        ax.set_title("Library depth", loc="left")
    else:
        _empty(ax, "no HiCRep contact-depth records")

    # (c) Pairtools duplicate candidates. Pairtools does not distinguish PCR
    # amplification from optical duplication, so the axis does not claim it does.
    ax = fig.add_subplot(gs[0, 2])
    _panel_label(ax, "c")
    d = lib.dropna(subset=["dup_pct"]).sort_values("dup_pct")
    if len(d):
        bars = ax.barh(range(len(d)), d["dup_pct"],
                       color=[_colour(c) for c in d["cell_type"]],
                       edgecolor="black", linewidth=0.3)
        failed = qc_gate_failures(d["dup_pct"], maximum=duplicate_gate)
        for bar, is_failed in zip(bars, failed):
            if is_failed:
                bar.set_hatch("///")
                bar.set_edgecolor("#b00020")
        ax.axvline(duplicate_gate, ls="--", lw=0.7, c="#b00020")
        _label_axis(ax, d.index, axis="y")
        ax.set_xlabel("Pairtools duplicate pairs (%)")
        ax.set_title("Library complexity", loc="left")
    else:
        _empty(ax, "no duplication stats")

    # (d) cis fraction.
    ax = fig.add_subplot(gs[1, 0])
    _panel_label(ax, "d")
    d = lib.dropna(subset=["cis_pct"]).sort_values("cis_pct", ascending=False)
    if len(d):
        bars = ax.bar(range(len(d)), d["cis_pct"],
                      color=[_colour(c) for c in d["cell_type"]],
                      edgecolor="black", linewidth=0.3)
        failed = qc_gate_failures(d["cis_pct"], minimum=cis_gate)
        for bar, is_failed in zip(bars, failed):
            if is_failed:
                bar.set_hatch("///")
                bar.set_edgecolor("#b00020")
        ax.axhline(cis_gate, ls="--", lw=0.7, c="#b00020")
        ax.set_ylim(0, 100)
        _label_axis(ax, d.index, axis="x")
        ax.set_ylabel("cis contacts (%)")
        ax.set_title("Cis / trans ratio", loc="left")
    else:
        _empty(ax, "no cis/trans stats")

    # (e) P(s) — contact frequency should generally decay with genomic distance.
    ax = fig.add_subplot(gs[1, 1:])
    _panel_label(ax, "e")
    drawn = 0
    for sid, r in lib.iterrows():
        f = results / f"qc/expected/{sid}.expected.cis.tsv"
        e = _required_tsv(
            f,
            f"expected-cis table for {sid}",
            {"dist_bp", "count.avg"},
        )
        # Use one normalization for every curve. A balance-aware expected-cis
        # rule may emit balanced values for one library and a raw fallback for
        # another; mixing those columns on one axis makes vertical differences
        # uninterpretable. Raw count averages are always the common denominator.
        col = "count.avg"
        # dist_bp, NOT dist. cooltools emits both: `dist` is a count of BINS and
        # `dist_bp` is that count times the bin size. Plotting `dist` under an axis
        # labelled "genomic separation (bp)" understated every separation by the bin
        # size -- at 25 kb the curve ran to 1e4 "bp" when it actually spans ~250 Mb.
        # The shape of P(s) is unchanged by the rescaling, which is exactly why this
        # was invisible: the figure looked right.
        dist_col = "dist_bp"
        g = e.groupby(dist_col)[col].mean()
        # Very long-distance bins contain few chromosome-region observations and
        # become visibly jagged. The QC question is whether distance decay behaves
        # sensibly over the informative HiChIP range, not whether a handful
        # of 100-Mb pixels happen to be non-zero.
        g = g[(g.index >= 25_000) & (g.index <= 50_000_000) & (g > 0)]
        if g.empty:
            continue
        ax.plot(g.index, g.values, lw=0.8, color=_colour(r["cell_type"]),
                alpha=0.85, ls="-" if r["hicrep_depth_adequate"] else ":")
        drawn += 1
    if drawn:
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("genomic separation (bp)")
        ax.set_ylabel("raw mean contact count")
        ax.set_title("Raw distance decay P(s)", loc="left")
        handles = [Line2D([0], [0], color=c, lw=1.2, label=k) for k, c in PALETTE.items()]
        handles.append(Line2D([0], [0], color=GREY, lw=1.2, ls=":", label="below depth floor"))
        ax.legend(handles=handles, frameon=False, loc="lower left")
    else:
        _empty(ax, "no expected-cis tables")

    fig.suptitle("Figure 1 · HiChIP library and contact-map quality", x=0.02,
                 ha="left", fontsize=9, fontweight="bold")
    # Colour is the only thing distinguishing cell types in panels a-d, so the key
    # belongs on the figure rather than inside one panel that may render empty.
    key = [Line2D([0], [0], color=c, lw=4, label=k)
           for k, c in PALETTE.items() if k in set(lib["cell_type"])]
    key.append(Line2D([0], [0], color="white", markerfacecolor="white",
                      markeredgecolor="#b00020", marker="s", markersize=7,
                      lw=0, label="below applicable configured gate (hatched)"))
    fig.legend(handles=key, frameon=False, ncol=len(key), fontsize=6,
               loc="upper right", bbox_to_anchor=(0.99, 1.0))
    _save(fig, out)


# ---------------------------------------------------------------- figure 2
def figure2_reproducibility(
    lib: pd.DataFrame, results: Path, out: Path, hicrep_threshold: float
) -> None:
    """Exploratory within-mark E1 similarity plus depth-aware HiCRep QC."""
    fig = plt.figure(figsize=(7.6, 3.6))
    # wspace: the colourbar label of panel a is drawn to its right and lands on
    # panel b's bar labels at anything tighter.
    gs = fig.add_gridspec(1, 2, wspace=0.95, width_ratios=[1.25, 1])

    # (a) E1 similarity is shown only within the same immunoprecipitated mark.
    # Cross-mark HiChIP E1 is not replicate QC or an unbiased compartment test.
    ax = fig.add_subplot(gs[0, 0])
    _panel_label(ax, "a")
    eigs = {}
    for sid in lib.index:
        f = results / f"qc/compartments/{sid}.cis.eigs.tsv"
        e = _required_tsv(
            f, f"compartment table for {sid}", {"chrom", "start", "E1"}
        )
        e = e.dropna(subset=["E1"])
        key = e["chrom"].astype(str) + ":" + e["start"].astype(str)
        eigs[sid] = pd.Series(e["E1"].values, index=key.values)

    if len(eigs) >= 2:
        # Pairwise-complete correlation preserves bins shared by each sample pair.
        # A global dropna() would keep only bins present in every library and can
        # make one sparse sample change every other pair's estimate.
        E = pd.DataFrame(eigs)
        C = E.corr(method="pearson", min_periods=100)
        same_mark = np.equal.outer(
            lib.loc[C.index, "mark"].astype(str).to_numpy(),
            lib.loc[C.columns, "mark"].astype(str).to_numpy(),
        )
        values = C.to_numpy(copy=True)
        values[~same_mark] = np.nan
        cmap = plt.get_cmap("RdBu_r").copy()
        cmap.set_bad("#eeeeee")
        im = ax.imshow(np.ma.masked_invalid(values), cmap=cmap, vmin=-1, vmax=1)
        _label_axis(ax, C.columns, axis="x")
        _label_axis(ax, C.index, axis="y")
        for i in range(len(C)):
            for j in range(len(C)):
                if np.isfinite(values[i, j]):
                    ax.text(j, i, f"{values[i, j]:.2f}", ha="center", va="center",
                            fontsize=4.5, color="black" if abs(values[i, j]) < 0.6 else "white")
        fig.colorbar(im, ax=ax, fraction=0.046, label="Pearson r (E1)")
        ax.set_title("Within-mark exploratory E1 similarity (100 kb)", loc="left")
    else:
        _empty(ax, "fewer than two compartment tables")

    # (b) HiCRep SCC, with depth-confounded pairs called out
    ax = fig.add_subplot(gs[0, 1])
    _panel_label(ax, "b")
    pairs = []
    seen = set()
    for sid in lib.index:
        f = results / f"qc/hicrep/{sid}.hicrep.json"
        h = _required_json(
            f, f"HiCRep QC for {sid}", {"status", "pairwise_scc", "contacts"}
        )
        for p in h.get("pairwise_scc", []):
            missing = sorted({"a", "b", "scc"} - set(p))
            if missing:
                raise ValueError(
                    f"HiCRep pair record for {sid} lacks fields {missing}: {f}"
                )
            key = tuple(sorted([p["a"], p["b"]]))
            if key in seen:
                continue
            seen.add(key)
            pairs.append(p)
    if pairs:
        pairs = sorted(pairs, key=lambda p: p["scc"])
        y = range(len(pairs))
        cols = ["#b00020" if p.get("depth_confounded") else "#4C72B0" for p in pairs]
        ax.barh(list(y), [p["scc"] for p in pairs], color=cols,
                edgecolor="black", linewidth=0.3)
        ax.set_yticks(list(y))
        ax.set_yticklabels([f"{p['a'].split('.')[0]}\nvs {p['b'].split('.')[0]}" for p in pairs],
                           fontsize=4.5)
        ax.set_xlabel("HiCRep SCC (25 kb)")
        # SCC can be negative. Clipping the axis at zero hides the strongest
        # possible evidence of irreproducibility.
        ax.set_xlim(-1, 1)
        ax.axvline(0, color="#666666", lw=0.5)
        ax.axvline(hicrep_threshold, color="#b00020", ls="--", lw=0.7)
        ax.set_title("Replicate concordance", loc="left")
        handles = [
            Line2D([0], [0], color="#4C72B0", lw=4, label="depth-adequate"),
            Line2D([0], [0], color="#b00020", lw=4, label="depth-confounded"),
            Line2D([0], [0], color="#b00020", lw=0.7, ls="--",
                   label=f"SCC gate {hicrep_threshold:g}"),
        ]
        # Keep the key outside the data rectangle: SCC bars can legitimately span
        # the full -1..1 range, so no in-panel corner is guaranteed to be empty.
        ax.legend(
            handles=handles,
            frameon=False,
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0,
        )
    else:
        _empty(ax, "no HiCRep results")

    fig.suptitle("Figure 2 · Reproducibility", x=0.02, ha="left",
                 fontsize=9, fontweight="bold")
    _save(fig, out)


# ---------------------------------------------------------------- figure 3
def figure3_loops(
    lib: pd.DataFrame, results: Path, out: Path, apa_bin_size: int
) -> None:
    primary = primary_reporting_libraries(lib)
    demonstration = lib[lib["report_role"] == "demonstration"].copy()
    if primary.empty:
        raise ValueError("Figure 3 requires at least one primary reporting library")
    fig = plt.figure(figsize=(7.6, 5.4))
    gs = fig.add_gridspec(2, 3, hspace=1.0, wspace=0.45)

    # (a) loop counts
    ax = fig.add_subplot(gs[0, 0])
    _panel_label(ax, "a")
    d = primary.sort_values("n_loops", ascending=False)
    if d["n_loops"].sum() > 0:
        ax.bar(range(len(d)), d["n_loops"],
               color=[_colour(c) for c in d["cell_type"]],
               edgecolor="black", linewidth=0.3)
        _label_axis(ax, d.index, axis="x")
        ax.set_ylabel("significant loops")
        ax.set_title("FitHiChIP loops", loc="left")
    else:
        _empty(ax, "no loops called")

    # (b) loop yield against contacts in the actual FitHiChIP distance range.
    # Marker area exposes anchor opportunity; total UU pairs would mix trans,
    # diagonal, and out-of-range contacts that the caller cannot turn into loops.
    ax = fig.add_subplot(gs[0, 1])
    _panel_label(ax, "b")
    d = primary.dropna(subset=["loop_search_contacts"])
    if len(d) and d["n_loops"].sum() > 0:
        anchor_counts = pd.to_numeric(d["n_consensus_peaks"], errors="coerce")
        anchor_scale = max(float(anchor_counts.max()), 1.0)
        for ct, sub in d.groupby("cell_type"):
            sizes = 14 + 36 * (
                pd.to_numeric(sub["n_consensus_peaks"], errors="coerce")
                .fillna(0) / anchor_scale
            )
            ax.scatter(sub["loop_search_contacts"] / 1e6, sub["n_loops"], s=sizes,
                       color=_colour(ct), edgecolor="black", linewidth=0.3, label=ct)
        ax.set_xscale("log")
        ax.set_xlabel("cis contacts in FitHiChIP range (millions)")
        ax.set_ylabel("significant loops")
        ax.legend(frameon=False)
        ax.set_title("Loop yield vs depth", loc="left")
        ax.text(
            0.5, -0.35, "marker area = assay-stratum consensus-anchor opportunity",
            transform=ax.transAxes, fontsize=5, color=GREY, ha="center", va="top",
        )
    else:
        _empty(ax, "no loops called")

    # (c) loop span distribution
    ax = fig.add_subplot(gs[0, 2])
    _panel_label(ax, "c")
    drawn = 0
    for sid, r in primary.iterrows():
        f = results / f"loops/{sid}/{sid}.interactions_FitHiC_{Q_LABEL}.bed"
        lp = _required_pipeline_loops(f, f"FitHiChIP loops for {sid}")
        if lp.empty:
            continue
        span = (lp["start2"].astype(float) - lp["start1"].astype(float)).abs()
        span = span[span > 0]
        if span.empty:
            continue
        ax.hist(np.log10(span), bins=40, histtype="step", lw=0.8,
                color=_colour(r["cell_type"]), density=True)
        drawn += 1
    if drawn:
        ax.set_xlabel("log10 loop span (bp)")
        ax.set_ylabel("density")
        ax.set_title("Loop span", loc="left")
    else:
        _empty(ax, "no loops called")

    # (d) APA, aggregated per (cell type, mark)
    by_ct: dict[tuple[str, str], list[np.ndarray]] = {}
    raw_fallbacks: dict[tuple[str, str], int] = {}
    for sid, r in primary.iterrows():
        f = results / f"qc/apa/{sid}.apa.npy"
        meta_f = results / f"qc/apa/{sid}.apa.json"
        _required_file(f, f"APA matrix for {sid}")
        apa_meta = _required_json(
            meta_f,
            f"APA metadata for {sid}",
            {"status", "normalization", "balance_converged"},
        )
        try:
            m = np.load(f, allow_pickle=False)
        except Exception as exc:
            raise ValueError(f"APA matrix is unreadable for {sid}: {f}") from exc
        status = validate_apa_matrix_contract(apa_meta, m, str(sid))
        if status == "NOT_ASSESSED":
            continue

        # Keyed on (cell_type, mark), NOT cell type alone. A TF anchor set and a
        # histone-mark anchor set have different APA profiles by construction.
        key = (r["cell_type"], r["mark"])
        # Raw fallbacks remain useful per-sample diagnostics, but their scale is
        # not commensurate with balanced matrices. Never average the two.
        if (
            apa_meta.get("balance_converged") is True
            and apa_meta.get("normalization") == "ICE-balanced"
        ):
            by_ct.setdefault(key, []).append(m)
        else:
            raw_fallbacks[key] = raw_fallbacks.get(key, 0) + 1

    # Preserve sample-sheet order while keeping each biological condition/mark
    # stratum separate. This remains correct when users replace the bundled cohort.
    panels = list(
        primary.reset_index()[["cell_type", "mark"]]
        .drop_duplicates()
        .itertuples(index=False, name=None)
    )
    fig.set_size_inches(max(7.6, 1.8 * len(panels)), 5.4)
    sub = gs[1, :].subgridspec(1, len(panels), wspace=0.62)

    for i, key in enumerate(panels):
        ct, mark = key
        ax = fig.add_subplot(sub[0, i])
        if i == 0:
            _panel_label(ax, "d")
        mats = by_ct.get(key)
        if not mats:
            excluded = raw_fallbacks.get(key, 0)
            if excluded:
                _empty(
                    ax,
                    f"{ct} {mark}: no balanced APA\n"
                    f"({excluded} raw-count fallback excluded)",
                )
            else:
                _empty(ax, f"{ct} {mark}: no assessed APA")
            continue
        m = np.mean(mats, axis=0)
        win = (m.shape[0] - 1) // 2
        im = ax.imshow(np.log2(m + 1), cmap="Reds", origin="lower",
                       extent=[-win, win, -win, win])
        excluded = raw_fallbacks.get(key, 0)
        suffix = f"; {excluded} raw excluded" if excluded else ""
        ax.set_title(f"{ct} {mark} (n={len(mats)} balanced{suffix})", loc="left",
                     color=_colour(ct), fontsize=7)
        bin_label = f"bins ({apa_bin_size / 1000:g} kb)"
        ax.set_xlabel(bin_label)
        if i == 0:
            ax.set_ylabel(bin_label)
        else:
            ax.set_yticklabels([])
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.ax.tick_params(labelsize=5)
        # Only the last colourbar is labelled: the label is drawn to the right of the
        # bar and ran straight into the next panel's tick labels.
        if i == len(panels) - 1:
            cb.set_label("log2(1+contacts)", fontsize=6)

    fig.suptitle("Figure 3 · Loop calling and aggregate peak analysis", x=0.02,
                 ha="left", fontsize=9, fontweight="bold")
    if len(demonstration):
        fig.text(
            0.5,
            0.01,
            f"{len(demonstration)} demonstration libraries are excluded from all "
            "Figure 3 headline panels; they remain in QC tables.",
            ha="center",
            va="bottom",
            fontsize=5.5,
            color=GREY,
        )
    _save(fig, out)


# ---------------------------------------------------------------- figure 4
DIFFERENTIAL_COLUMNS = {
    "loop_key", "chrom1", "start1", "end1", "chrom2", "start2", "end2", "padj",
    "log2FoldChange", "lfcSE", "analysis_status",
}
PAIRED_EFFECT_COLUMNS = {
    "loop_key", "chrom1", "start1", "end1", "chrom2", "start2", "end2",
    "pair_id", "case_normalized_count", "control_normalized_count",
    "paired_log2_ratio", "analysis_status",
}
DIFFERENTIAL_DESIGN_KEYS = {
    "comparison", "analysis_status", "n_complete_pairs",
    "publication_eligible", "publication_min_complete_pairs", "paired_subjects",
    "candidate_loops", "tested_loops",
}


def load_differential_bundle(
    comparison: str, results: Path
) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    """Load and cross-check the three required Figure 4 contracts."""
    base = results / f"diff/{comparison}"
    result = _required_tsv(
        base / "differential_loops.tsv",
        f"differential results for {comparison}",
        DIFFERENTIAL_COLUMNS,
        allow_empty_rows=False,
    )
    design = _required_json(
        base / "design.json",
        f"differential design for {comparison}",
        DIFFERENTIAL_DESIGN_KEYS,
    )
    paired = _required_tsv(
        base / "paired_effects.tsv",
        f"paired effects for {comparison}",
        PAIRED_EFFECT_COLUMNS,
        allow_empty_rows=False,
    )
    if str(design["comparison"]) != comparison:
        raise ValueError(
            f"differential design comparison mismatch: expected {comparison!r}, "
            f"got {design['comparison']!r}"
        )
    status = str(design["analysis_status"])
    allowed = {"PILOT_UNDERPOWERED", "STANDARD_INFERENCE"}
    if status not in allowed:
        raise ValueError(
            f"differential analysis_status must be one of {sorted(allowed)}; got {status!r}"
        )
    if type(design["publication_eligible"]) is not bool:
        raise ValueError("publication_eligible must be a JSON boolean")
    n_complete_pairs = design["n_complete_pairs"]
    if type(n_complete_pairs) is not int or n_complete_pairs < 1:
        raise ValueError("n_complete_pairs must be a positive JSON integer")
    publication_min = design["publication_min_complete_pairs"]
    if type(publication_min) is not int or publication_min < 2:
        raise ValueError(
            "publication_min_complete_pairs must be a JSON integer of at least two"
        )
    tested_loops = design["tested_loops"]
    candidate_loops = design["candidate_loops"]
    if type(tested_loops) is not int or tested_loops < 1:
        raise ValueError("tested_loops must be a positive JSON integer")
    if type(candidate_loops) is not int or candidate_loops < tested_loops:
        raise ValueError("candidate_loops must be an integer >= tested_loops")
    if len(result) != tested_loops:
        raise ValueError(
            f"differential result has {len(result)} rows but design declares "
            f"tested_loops={tested_loops}"
        )
    subjects = design["paired_subjects"]
    if not isinstance(subjects, list) or any(
        not isinstance(subject, str) or not subject for subject in subjects
    ):
        raise ValueError("paired_subjects must be a list of non-empty strings")
    expected_pairs = set(subjects)
    if len(expected_pairs) != n_complete_pairs or len(subjects) != n_complete_pairs:
        raise ValueError("paired_subjects does not match n_complete_pairs")

    for label, frame in (("differential results", result), ("paired effects", paired)):
        if frame["analysis_status"].isna().any():
            raise ValueError(f"{label} contains null analysis_status values")
        observed = set(frame["analysis_status"].dropna().astype(str).unique())
        if observed and observed != {status}:
            raise ValueError(
                f"{label} status {sorted(observed)} disagrees with design status {status}"
            )
    if design["publication_eligible"] != (status == "STANDARD_INFERENCE"):
        raise ValueError(
            "publication_eligible disagrees with differential analysis_status"
        )
    if result["loop_key"].astype(str).duplicated().any():
        raise ValueError("differential results contain duplicate loop_key rows")
    if result["loop_key"].isna().any() or paired[["loop_key", "pair_id"]].isna().any().any():
        raise ValueError("differential loop_key and pair_id values must be non-null")
    result = _validate_loop_coordinates(
        result, f"differential results for {comparison}"
    )
    numeric_contracts = {
        "padj": (0.0, 1.0),
        "log2FoldChange": (None, None),
        "lfcSE": (0.0, None),
    }
    for column, (lower, upper) in numeric_contracts.items():
        try:
            values = pd.to_numeric(result[column], errors="raise")
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"differential results have nonnumeric {column} values"
            ) from exc
        present = values.notna()
        finite = values.loc[present].to_numpy(dtype=float)
        if not np.isfinite(finite).all():
            raise ValueError(
                f"differential results have non-finite {column} values"
            )
        if lower is not None and (values.loc[present] < lower).any():
            raise ValueError(
                f"differential results have {column} below {lower:g}"
            )
        if upper is not None and (values.loc[present] > upper).any():
            raise ValueError(
                f"differential results have {column} above {upper:g}"
            )
        result[column] = values.astype(float)
    paired = _validate_loop_coordinates(
        paired, f"paired effects for {comparison}"
    )
    paired_keys = paired[["loop_key", "pair_id"]].astype(str)
    if paired_keys.duplicated().any():
        raise ValueError("paired effects contain duplicate (loop_key, pair_id) rows")
    result_loops = set(result["loop_key"].astype(str))
    paired_loops = set(paired["loop_key"].astype(str))
    if paired_loops != result_loops:
        raise ValueError("paired-effect loop coverage does not match differential results")
    coordinate_columns = [
        "chrom1", "start1", "end1", "chrom2", "start2", "end2"
    ]
    paired_coordinates = paired[
        ["loop_key", *coordinate_columns]
    ].drop_duplicates()
    if len(paired_coordinates) != tested_loops:
        raise ValueError("paired-effect coordinates are inconsistent within a loop_key")
    expected_coordinates = result.set_index("loop_key")[coordinate_columns].sort_index()
    observed_coordinates = paired_coordinates.set_index("loop_key")[
        coordinate_columns
    ].sort_index()
    if not expected_coordinates.equals(observed_coordinates):
        raise ValueError("paired-effect coordinates disagree with differential results")
    expected_paired_rows = tested_loops * n_complete_pairs
    if len(paired) != expected_paired_rows:
        raise ValueError(
            f"paired effects has {len(paired)} rows but expected "
            f"{tested_loops} x {n_complete_pairs} = {expected_paired_rows}"
        )
    observed_pairs = set(paired["pair_id"].astype(str))
    if observed_pairs != expected_pairs:
        raise ValueError("paired-effect pair_id values do not match paired_subjects")
    coverage = paired.assign(
        loop_key=paired["loop_key"].astype(str),
        pair_id=paired["pair_id"].astype(str),
    ).groupby("loop_key")["pair_id"].agg(set)
    if any(value != expected_pairs for value in coverage):
        raise ValueError("one or more loops lack complete paired-subject coverage")
    for column in (
        "case_normalized_count", "control_normalized_count", "paired_log2_ratio"
    ):
        try:
            values = pd.to_numeric(paired[column], errors="raise")
        except (TypeError, ValueError) as exc:
            raise ValueError(f"paired effects has nonnumeric {column}") from exc
        if values.isna().any() or not np.isfinite(values.to_numpy(dtype=float)).all():
            raise ValueError(f"paired effects has non-finite {column}")
        paired[column] = values.astype(float)
    if (paired[["case_normalized_count", "control_normalized_count"]] < 0).any().any():
        raise ValueError("paired normalized counts must be non-negative")
    return result, design, paired


def _forest_label(row: pd.Series) -> str:
    return (
        f"{row['chrom1']}:{int(row['start1']):,}\n"
        f"{row['chrom2']}:{int(row['start2']):,}"
    )


def figure4_differential(comparisons: list[str], results: Path, out: Path,
                         fdr: float, lfc: float) -> None:
    n = max(len(comparisons), 1)
    fig = plt.figure(figsize=(max(3.0 * n, 5.2), 8.2))
    gs = fig.add_gridspec(3, n, hspace=0.72, wspace=0.55)

    for i, comp in enumerate(comparisons):
        res, design, paired = load_differential_bundle(comp, results)
        status = str(design["analysis_status"])
        pilot = status == "PILOT_UNDERPOWERED"
        status_line = (
            f"{status}; {int(design['n_complete_pairs'])} complete pairs"
        )

        # (a) Model-wide overview. Pilot results retain their model values but are
        # explicitly labelled exploratory rather than publication-ready hits.
        ax = fig.add_subplot(gs[0, i])
        if i == 0:
            _panel_label(ax, "a")
        finite = res.dropna(subset=["padj", "log2FoldChange"]).copy()
        if finite.empty:
            _empty(ax, f"{comp}\nno finite model estimates\n{status_line}")
        else:
            sig = (
                (finite["padj"] < fdr)
                & (finite["log2FoldChange"].abs() >= lfc)
            )
            y = -np.log10(finite["padj"].clip(lower=np.nextafter(0, 1)))
            ax.scatter(
                finite.loc[~sig, "log2FoldChange"], y[~sig], s=2,
                c="#d9d9d9", rasterized=True, linewidths=0,
            )
            up = sig & (finite["log2FoldChange"] > 0)
            down = sig & (finite["log2FoldChange"] < 0)
            qualifier = "exploratory " if pilot else ""
            ax.scatter(
                finite.loc[up, "log2FoldChange"], y[up], s=3, c="#C44E52",
                rasterized=True, linewidths=0,
                label=f"{qualifier}up ({int(up.sum())})",
            )
            ax.scatter(
                finite.loc[down, "log2FoldChange"], y[down], s=3, c="#4C72B0",
                rasterized=True, linewidths=0,
                label=f"{qualifier}down ({int(down.sum())})",
            )
            ax.axhline(-np.log10(fdr), ls="--", lw=0.5, c="black")
            ax.axvline(lfc, ls="--", lw=0.5, c="black")
            ax.axvline(-lfc, ls="--", lw=0.5, c="black")
            ax.set_xlabel("model log2 fold change")
            if i == 0:
                ax.set_ylabel("-log10 adjusted p")
            ax.legend(frameon=False, loc="upper left")
        ax.set_title(
            f"{comp.replace('_', ' ')}\n{status_line}",
            loc="left",
            color="#b00020" if pilot else "black",
        )

        # (b) Donor-level paired effects for the best-ranked tested loops.
        ax2 = fig.add_subplot(gs[1, i])
        if i == 0:
            _panel_label(ax2, "b")
        ranked = finite.sort_values("padj", kind="stable").head(12)
        top_keys = ranked["loop_key"].astype(str).tolist()
        donor = paired[paired["loop_key"].astype(str).isin(top_keys)].copy()
        if top_keys and not donor.empty:
            heat = donor.pivot(
                index="loop_key", columns="pair_id", values="paired_log2_ratio"
            ).reindex(top_keys)
            heat = heat.dropna(how="all")
        else:
            heat = pd.DataFrame()
        if not heat.empty:
            limit = float(np.nanmax(np.abs(heat.to_numpy(dtype=float))))
            limit = max(limit, 0.5)
            im = ax2.imshow(
                heat.to_numpy(dtype=float), aspect="auto", cmap="RdBu_r",
                vmin=-limit, vmax=limit,
            )
            ax2.set_xticks(range(len(heat.columns)))
            ax2.set_xticklabels(heat.columns, rotation=45, ha="right", fontsize=5)
            ax2.set_yticks(range(len(heat.index)))
            ax2.set_yticklabels([str(key)[:8] for key in heat.index], fontsize=4.5)
            ax2.set_xlabel("paired donor / subject")
            if i == 0:
                ax2.set_ylabel("top loop key")
            fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.04,
                         label="paired log2 ratio")
            ax2.set_title("Within-pair normalized effects", loc="left")
        else:
            _empty(ax2, "no paired effects for ranked loops")

        # (c) Model uncertainty: effect +/- 1.96 * lfcSE. These are Wald
        # intervals, not bootstrap donor intervals, and pilot panels say so.
        ax3 = fig.add_subplot(gs[2, i])
        if i == 0:
            _panel_label(ax3, "c")
        forest = res.dropna(subset=["padj", "log2FoldChange", "lfcSE"]).copy()
        forest = forest[forest["lfcSE"] >= 0].sort_values("padj", kind="stable").head(10)
        if forest.empty:
            _empty(ax3, "no finite Wald intervals")
        else:
            forest = forest.iloc[::-1].reset_index(drop=True)
            y_pos = np.arange(len(forest))
            effect = forest["log2FoldChange"].to_numpy(dtype=float)
            interval = 1.96 * forest["lfcSE"].to_numpy(dtype=float)
            significant = (
                (forest["padj"].to_numpy(dtype=float) < fdr)
                & (np.abs(effect) >= lfc)
            )
            for j, (estimate, error, is_sig) in enumerate(
                zip(effect, interval, significant)
            ):
                ax3.errorbar(
                    estimate, j, xerr=error, fmt="o", markersize=2.8,
                    color="#C44E52" if is_sig else "#666666",
                    ecolor="#C44E52" if is_sig else "#999999",
                    elinewidth=0.7, capsize=1.5,
                )
            ax3.axvline(0, color="black", lw=0.6)
            ax3.set_yticks(y_pos)
            ax3.set_yticklabels(
                [_forest_label(row) for _, row in forest.iterrows()], fontsize=4.5
            )
            ax3.set_xlabel("log2 fold change (95% Wald interval)")
            ax3.set_title(
                "Exploratory model intervals" if pilot else "Model effect intervals",
                loc="left", color="#b00020" if pilot else "black",
            )

    if not comparisons:
        for row, letter in enumerate(("a", "b", "c")):
            ax = fig.add_subplot(gs[row, 0])
            _panel_label(ax, letter)
            _empty(ax, "no differential comparisons configured")

    fig.suptitle("Figure 4 · Paired differential loops", x=0.02, ha="left",
                 fontsize=9, fontweight="bold")
    _save(fig, out)


# ---------------------------------------------------------------- figure 5
STRIPE_COLUMNS = {
    "chr", "pos1", "pos2", "chr2", "pos3", "pos4", "length", "width",
    "Mean", "maxpixel", "pvalue", "Stripiness",
}


def _required_stripe_table(path: str | Path, label: str) -> pd.DataFrame:
    frame = _required_tsv(path, label, STRIPE_COLUMNS)
    if frame.empty:
        return frame
    try:
        length = pd.to_numeric(frame["length"], errors="raise")
    except (TypeError, ValueError) as exc:
        raise ValueError(f"required {label} has nonnumeric stripe length") from exc
    if (
        length.isna().any()
        or not np.isfinite(length.to_numpy(dtype=float)).all()
        or (length <= 0).any()
    ):
        raise ValueError(f"required {label} has non-positive or non-finite length")
    frame["length"] = length.astype(float)
    return frame


def figure5_stripes(
    lib: pd.DataFrame, stripe_files: list[str | Path], out: Path
) -> None:
    """Stripes, read against the anchor they were called on.

    Only libraries configured with the ``primary`` report role enter this headline
    summary. Demonstration libraries remain in every upstream result, the library
    table, and the QC figures, but cannot invite an unsupported cross-mark claim.
    """
    fig = plt.figure(figsize=(7.6, 4.0))
    # top=0.82: the suptitle sits on the same line as the panel "a" label otherwise.
    gs = fig.add_gridspec(1, 3, wspace=0.5, top=0.82, bottom=0.32)

    primary = primary_reporting_libraries(lib)
    demonstration = lib[lib["report_role"] == "demonstration"]
    stripes: dict[str, pd.DataFrame] = {}
    for file in stripe_files:
        f = Path(file)
        df = _required_stripe_table(f, f"stripe table {f}")
        sid = f.parent.name
        if sid not in primary.index:
            continue
        stripes[sid] = df

    # (a) count by sample, coloured by mark
    ax = fig.add_subplot(gs[0, 0])
    _panel_label(ax, "a")
    if stripes:
        counts = pd.Series({s: len(d) for s, d in stripes.items()})
        counts = counts.reindex(primary.index).fillna(0).sort_values(ascending=False)
        cols = ["#8172B2" if primary.loc[s, "mark"] == "CTCF" else "#CCB974" for s in counts.index]
        ax.bar(range(len(counts)), counts.values, color=cols,
               edgecolor="black", linewidth=0.3)
        _label_axis(ax, counts.index, axis="x")
        ax.set_ylabel("stripes")
        ax.set_title("Stripe count", loc="left")
        marks = sorted({str(primary.loc[s, "mark"]) for s in counts.index})
        ax.legend(
            handles=[
                Line2D(
                    [0], [0],
                    color="#8172B2" if mark == "CTCF" else "#CCB974",
                    lw=4,
                    label=f"{mark} anchors",
                )
                for mark in marks
            ],
            frameon=False,
        )
    else:
        _empty(ax, "no stripe tables")

    # (b) stripe yield against depth, NOT stripes-per-mark.
    #
    # A CTCF-vs-H3K27ac dot plot invites an anchor-class interpretation, but anchor
    # class is easily confounded with depth and library complexity, while stripe count
    # tracks depth just as loop count does. The honest panel exposes that relationship
    # directly rather than reducing each mark to one mean.
    ax = fig.add_subplot(gs[0, 1])
    _panel_label(ax, "b")
    rec = [
        {
            "mark": primary.loc[s, "mark"], "n": len(d),
            "depth": primary.loc[s, "stripe_search_contacts"],
        }
        for s, d in stripes.items()
        if s in primary.index and not pd.isna(primary.loc[s, "stripe_search_contacts"])
    ]
    df = pd.DataFrame(rec)
    if not df.empty:
        for mark, sub in df.groupby("mark"):
            col = "#8172B2" if mark == "CTCF" else "#CCB974"
            ax.scatter(sub["depth"] / 1e6, sub["n"], s=22, color=col,
                       edgecolor="black", linewidth=0.3, label=f"{mark} anchors")
        ax.set_xscale("log")
        ax.set_xlabel("primary cis off-diagonal contacts (millions)")
        ax.set_ylabel("stripes per library")
        ax.set_title("Stripe yield vs depth", loc="left")
        ax.legend(frameon=False, loc="upper left")
        ax.text(0.5, -0.42,
                "Compare marks only after checking depth and complexity;\n"
                "this panel makes potential confounding visible.",
                transform=ax.transAxes, fontsize=5, color=GREY, ha="center", va="top")
    else:
        _empty(ax, "no stripe tables")

    # (c) stripe length
    ax = fig.add_subplot(gs[0, 2])
    _panel_label(ax, "c")
    drawn = 0
    for sid, d in stripes.items():
        if d.empty:
            continue
        L = d["length"]
        if len(L) < 2:
            continue
        col = "#8172B2" if primary.loc[sid, "mark"] == "CTCF" else "#CCB974"
        ax.hist(np.log10(L), bins=30, histtype="step", lw=0.8, color=col, density=True)
        drawn += 1
    if drawn:
        ax.set_xlabel("log10 stripe length (bp)")
        ax.set_ylabel("density")
        ax.set_title("Stripe length", loc="left")
    else:
        _empty(ax, "no stripes called")

    fig.suptitle("Figure 5 · Architectural stripes", x=0.02, ha="left",
                 fontsize=9, fontweight="bold")
    if len(demonstration):
        fig.text(
            0.5,
            0.03,
            f"{len(demonstration)} configured demonstration libraries are excluded "
            "from this headline summary; they remain visible in QC and raw tables.",
            ha="center",
            va="bottom",
            fontsize=5.5,
            color=GREY,
        )
    _save(fig, out)


# ---------------------------------------------------------------- main
def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    _require_plotting()

    # Snakemake guarantees dependency ordering, while this explicit check makes
    # the direct script contract equally clear and catches broken symlinks.
    for input_path in list(snakemake.input):
        if not Path(input_path).is_file():
            raise FileNotFoundError(f"required publication-figure input is missing: {input_path}")

    global Q_LABEL
    Q_LABEL = snakemake.params.q_label

    results = Path(snakemake.params.results)
    outdir = Path(snakemake.params.outdir)
    # A clean first run has no results/figures directory yet. Create it before
    # writing the cohort table; the plotting helper only creates it later.
    outdir.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(
        snakemake.input.samples,
        sep="\t",
        comment="#",
        dtype=str,
        keep_default_na=False,
    )
    required_sample_columns = {"sample_id", "cell_type", "mark"}
    missing_sample_columns = sorted(required_sample_columns - set(samples.columns))
    if missing_sample_columns:
        raise ValueError(
            f"samples TSV lacks figure metadata columns: {missing_sample_columns}"
        )
    if samples["sample_id"].duplicated().any():
        raise ValueError("samples TSV contains duplicate sample_id values")
    samples = samples.set_index("sample_id")
    min_contacts = int(snakemake.params.min_contacts)
    fdr = float(snakemake.params.differential_fdr)
    lfc = float(snakemake.params.differential_log2fc_min)
    comparisons = list(snakemake.params.comparisons)

    lib = _library_table(
        samples,
        results,
        min_contacts,
        set(snakemake.params.demonstration_samples),
    )
    lib.to_csv(snakemake.output.table, sep="\t")
    log.info("library table:\n%s", lib.to_string())

    figure1_library_qc(
        lib,
        results,
        outdir / "figure1_library_qc",
        min_contacts,
        dict(snakemake.params.qc_thresholds),
    )
    figure2_reproducibility(
        lib,
        results,
        outdir / "figure2_reproducibility",
        float(snakemake.params.hicrep_threshold),
    )
    figure3_loops(
        lib,
        results,
        outdir / "figure3_loops_apa",
        int(snakemake.params.apa_bin_size),
    )
    figure4_differential(comparisons, results, outdir / "figure4_differential", fdr, lfc)
    figure5_stripes(
        lib, list(snakemake.input.stripes), outdir / "figure5_stripes"
    )


# Guarded so the module can be imported by the tests. Snakemake injects
# `snakemake` into the script's globals before executing it.
if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
