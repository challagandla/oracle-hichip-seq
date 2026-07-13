"""
Cohort-level publication figures.

The per-sample PNGs the other stages emit are diagnostics. These are the figures
that carry the argument, so they are assembled across the whole cohort, written as
vector PDF alongside 400 dpi PNG, and every panel that can be empty renders an
explicit "no data" instead of taking the stage down -- a single loopless or failed
library must not destroy the figure set at the end of a multi-day run.

Nothing here silently hides a bad library. Shallow libraries stay on the plots and
are marked, because the reader's first question about a HiChIP cohort is whether
the contrast is driven by biology or by depth.
"""
import json
import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_loops_bedpe, setup_logging  # noqa: E402

log = logging.getLogger(__name__)

# ---------------------------------------------------------------- style
plt.rcParams.update({
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
})

PALETTE = {"Naive": "#4C72B0", "Th17": "#C44E52", "Treg": "#55A868"}
GREY = "#9e9e9e"


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


def _empty(ax, msg: str) -> None:
    ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=6,
            color=GREY, transform=ax.transAxes, wrap=True)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def _save(fig, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"))
    fig.savefig(stem.with_suffix(".png"), dpi=400)
    plt.close(fig)
    log.info("wrote %s.{pdf,png}", stem)


# ---------------------------------------------------------------- parsing
def _pairtools_stats(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    if not path.exists():
        return out
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            try:
                out[parts[0]] = float(parts[1])
            except ValueError:
                pass
    return out


def _library_table(samples: pd.DataFrame, results: Path, min_contacts: int) -> pd.DataFrame:
    rows = []
    for sid, meta in samples.iterrows():
        dedup = _pairtools_stats(results / f"qc/pairtools/{sid}.dedup.stats.txt")
        total = dedup.get("total", np.nan)
        nodups = dedup.get("total_nodups", np.nan)
        dups = dedup.get("total_dups", np.nan)
        cis = dedup.get("cis", np.nan)
        trans = dedup.get("trans", np.nan)

        loops_f = results / f"loops/{sid}/{sid}.interactions_FitHiC_{Q_LABEL}.bed"
        try:
            n_loops = len(load_loops_bedpe(loops_f))
        except Exception:
            n_loops = 0

        apa_f = results / f"qc/apa/{sid}.apa.json"
        apa = json.loads(apa_f.read_text()) if apa_f.exists() else {}

        rows.append({
            "sample": sid,
            "cell_type": meta.get("cell_type", "?"),
            "mark": meta.get("mark", "?"),
            "unique_pairs": nodups,
            "dup_pct": 100.0 * dups / total if total and not np.isnan(total) and total > 0 else np.nan,
            "cis_pct": 100.0 * cis / (cis + trans) if (cis + trans) > 0 else np.nan,
            "n_loops": n_loops,
            "apa": apa.get("apa_vs_random_shift"),
            "usable": bool(nodups and nodups >= min_contacts),
        })
    return pd.DataFrame(rows).set_index("sample")


# ---------------------------------------------------------------- figure 1
def figure1_library_qc(lib: pd.DataFrame, results: Path, out: Path, min_contacts: int) -> None:
    fig = plt.figure(figsize=(7.2, 6.4))
    # Generous hspace: the sample labels are rotated 90 degrees under panels a and c,
    # and at the default spacing they run into the title of the panel below.
    gs = fig.add_gridspec(2, 2, hspace=1.05, wspace=0.45)

    # (a) unique pairs, with the depth floor drawn on
    ax = fig.add_subplot(gs[0, 0]); _panel_label(ax, "a")
    d = lib.dropna(subset=["unique_pairs"]).sort_values("unique_pairs", ascending=False)
    if len(d):
        colours = [_colour(c) for c in d["cell_type"]]
        bars = ax.bar(range(len(d)), d["unique_pairs"] / 1e6, color=colours,
                      edgecolor="black", linewidth=0.3)
        for i, (_, r) in enumerate(d.iterrows()):
            if not r["usable"]:
                bars[i].set_hatch("///")
                bars[i].set_edgecolor("#b00020")
        ax.axhline(min_contacts / 1e6, ls="--", lw=0.7, c="#b00020")
        # Annotate above the axes, not inside them: the line sits low on a log axis
        # and the label lands on top of the tick labels.
        ax.text(0.99, 1.01, f"{min_contacts/1e6:.0f}M usability floor",
                transform=ax.transAxes, fontsize=5, c="#b00020", ha="right", va="bottom")
        _label_axis(ax, d.index, axis="x")
        ax.set_ylabel("unique pairs (millions)")
        ax.set_yscale("log")
        ax.set_title("Library depth", loc="left")
    else:
        _empty(ax, "no pairtools stats")

    # (b) duplication
    ax = fig.add_subplot(gs[0, 1]); _panel_label(ax, "b")
    d = lib.dropna(subset=["dup_pct"]).sort_values("dup_pct")
    if len(d):
        ax.barh(range(len(d)), d["dup_pct"],
                color=[_colour(c) for c in d["cell_type"]],
                edgecolor="black", linewidth=0.3)
        ax.axvline(50, ls="--", lw=0.7, c="#b00020")
        _label_axis(ax, d.index, axis="y")
        ax.set_xlabel("PCR duplicates (%)")
        ax.set_title("Library complexity", loc="left")
    else:
        _empty(ax, "no duplication stats")

    # (c) cis fraction — the single best HiChIP sanity metric
    ax = fig.add_subplot(gs[1, 0]); _panel_label(ax, "c")
    d = lib.dropna(subset=["cis_pct"]).sort_values("cis_pct", ascending=False)
    if len(d):
        ax.bar(range(len(d)), d["cis_pct"],
               color=[_colour(c) for c in d["cell_type"]],
               edgecolor="black", linewidth=0.3)
        ax.axhline(70, ls="--", lw=0.7, c="#b00020")
        ax.set_ylim(0, 100)
        _label_axis(ax, d.index, axis="x")
        ax.set_ylabel("cis contacts (%)")
        ax.set_title("Cis / trans ratio", loc="left")
    else:
        _empty(ax, "no cis/trans stats")

    # (d) P(s) — distance decay; a real contact map falls off ~ s^-1
    ax = fig.add_subplot(gs[1, 1]); _panel_label(ax, "d")
    drawn = 0
    for sid, r in lib.iterrows():
        f = results / f"qc/expected/{sid}.expected.cis.tsv"
        if not f.exists():
            continue
        try:
            e = pd.read_csv(f, sep="\t")
        except Exception:
            continue
        col = "balanced.avg" if "balanced.avg" in e.columns else (
            "count.avg" if "count.avg" in e.columns else None)
        # dist_bp, NOT dist. cooltools emits both: `dist` is a count of BINS and
        # `dist_bp` is that count times the bin size. Plotting `dist` under an axis
        # labelled "genomic separation (bp)" understated every separation by the bin
        # size -- at 25 kb the curve ran to 1e4 "bp" when it actually spans ~250 Mb.
        # The shape of P(s) is unchanged by the rescaling, which is exactly why this
        # was invisible: the figure looked right.
        dist_col = "dist_bp" if "dist_bp" in e.columns else None
        if col is None or dist_col is None:
            continue
        g = e.groupby(dist_col)[col].mean()
        g = g[(g.index > 0) & (g > 0)]
        if g.empty:
            continue
        ax.plot(g.index, g.values, lw=0.8, color=_colour(r["cell_type"]),
                alpha=0.85, ls="-" if r["usable"] else ":")
        drawn += 1
    if drawn:
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("genomic separation (bp)")
        ax.set_ylabel("mean contact frequency")
        ax.set_title("Distance decay P(s)", loc="left")
        handles = [Line2D([0], [0], color=c, lw=1.2, label=k) for k, c in PALETTE.items()]
        handles.append(Line2D([0], [0], color=GREY, lw=1.2, ls=":", label="below depth floor"))
        ax.legend(handles=handles, frameon=False, loc="lower left")
    else:
        _empty(ax, "no expected-cis tables")

    fig.suptitle("Figure 1 · HiChIP library and contact-map quality", x=0.02,
                 ha="left", fontsize=9, fontweight="bold")
    # Colour is the only thing distinguishing cell types in panels a-c, so the key
    # belongs on the figure rather than inside one panel that may render empty.
    key = [Line2D([0], [0], color=c, lw=4, label=k)
           for k, c in PALETTE.items() if k in set(lib["cell_type"])]
    key.append(Line2D([0], [0], color="white", markerfacecolor="white",
                      markeredgecolor="#b00020", marker="s", markersize=7,
                      lw=0, label="below depth floor (hatched)"))
    fig.legend(handles=key, frameon=False, ncol=len(key), fontsize=6,
               loc="upper right", bbox_to_anchor=(0.99, 1.0))
    _save(fig, out)


# ---------------------------------------------------------------- figure 2
def figure2_reproducibility(lib: pd.DataFrame, results: Path, out: Path) -> None:
    """Compartment agreement is the honest all-vs-all reproducibility view here.

    HiCRep SCC is reported, but it is computed only within replicate groups and is
    depth-dominated (see rules/06_loop_qc.smk), so it cannot be laid out as a
    cohort-wide matrix without inviting exactly the misreading it is prone to.
    Compartment E1 at 100 kb is robust to depth and separates cell types.
    """
    fig = plt.figure(figsize=(7.6, 3.6))
    gs = fig.add_gridspec(1, 2, wspace=0.5, width_ratios=[1.2, 1])

    # (a) E1 correlation across all libraries
    ax = fig.add_subplot(gs[0, 0]); _panel_label(ax, "a")
    eigs = {}
    for sid in lib.index:
        f = results / f"qc/compartments/{sid}.cis.eigs.tsv"
        if not f.exists():
            continue
        try:
            e = pd.read_csv(f, sep="\t")
        except Exception:
            continue
        if "E1" not in e.columns:
            continue
        e = e.dropna(subset=["E1"])
        key = e["chrom"].astype(str) + ":" + e["start"].astype(str) if "chrom" in e.columns else None
        if key is None:
            continue
        eigs[sid] = pd.Series(e["E1"].values, index=key.values)

    if len(eigs) >= 2:
        E = pd.DataFrame(eigs).dropna()
        C = E.corr(method="pearson")
        im = ax.imshow(C.values, cmap="RdBu_r", vmin=-1, vmax=1)
        _label_axis(ax, C.columns, axis="x")
        _label_axis(ax, C.index, axis="y")
        for i in range(len(C)):
            for j in range(len(C)):
                ax.text(j, i, f"{C.values[i, j]:.2f}", ha="center", va="center",
                        fontsize=4.5, color="black" if abs(C.values[i, j]) < 0.6 else "white")
        fig.colorbar(im, ax=ax, fraction=0.046, label="Pearson r (E1)")
        ax.set_title("A/B compartment agreement (100 kb)", loc="left")
    else:
        _empty(ax, "fewer than two compartment tables")

    # (b) HiCRep SCC, with depth-confounded pairs called out
    ax = fig.add_subplot(gs[0, 1]); _panel_label(ax, "b")
    pairs = []
    seen = set()
    for sid in lib.index:
        f = results / f"qc/hicrep/{sid}.hicrep.json"
        if not f.exists():
            continue
        try:
            h = json.loads(f.read_text())
        except Exception:
            continue
        for p in h.get("pairwise_scc", []):
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
        ax.set_xlim(0, 1)
        ax.set_title("Replicate concordance", loc="left")
        handles = [
            Line2D([0], [0], color="#4C72B0", lw=4, label="depth-adequate"),
            Line2D([0], [0], color="#b00020", lw=4, label="depth-confounded"),
        ]
        ax.legend(handles=handles, frameon=False, loc="lower right")
    else:
        _empty(ax, "no HiCRep results")

    fig.suptitle("Figure 2 · Reproducibility", x=0.02, ha="left",
                 fontsize=9, fontweight="bold")
    _save(fig, out)


# ---------------------------------------------------------------- figure 3
def figure3_loops(lib: pd.DataFrame, results: Path, out: Path) -> None:
    fig = plt.figure(figsize=(7.6, 5.4))
    gs = fig.add_gridspec(2, 3, hspace=1.0, wspace=0.45)

    # (a) loop counts
    ax = fig.add_subplot(gs[0, 0]); _panel_label(ax, "a")
    d = lib.sort_values("n_loops", ascending=False)
    if d["n_loops"].sum() > 0:
        ax.bar(range(len(d)), d["n_loops"],
               color=[_colour(c) for c in d["cell_type"]],
               edgecolor="black", linewidth=0.3)
        _label_axis(ax, d.index, axis="x")
        ax.set_ylabel("significant loops")
        ax.set_title("FitHiChIP loops", loc="left")
    else:
        _empty(ax, "no loops called")

    # (b) loops vs depth — the plot that says whether loop count is biology or depth
    ax = fig.add_subplot(gs[0, 1]); _panel_label(ax, "b")
    d = lib.dropna(subset=["unique_pairs"])
    if len(d) and d["n_loops"].sum() > 0:
        for ct, sub in d.groupby("cell_type"):
            ax.scatter(sub["unique_pairs"] / 1e6, sub["n_loops"], s=14,
                       color=_colour(ct), edgecolor="black", linewidth=0.3, label=ct)
        ax.set_xscale("log")
        ax.set_xlabel("unique pairs (millions)")
        ax.set_ylabel("significant loops")
        ax.legend(frameon=False)
        ax.set_title("Loop yield vs depth", loc="left")
    else:
        _empty(ax, "no loops called")

    # (c) loop span distribution
    ax = fig.add_subplot(gs[0, 2]); _panel_label(ax, "c")
    drawn = 0
    for sid, r in lib.iterrows():
        f = results / f"loops/{sid}/{sid}.interactions_FitHiC_{Q_LABEL}.bed"
        try:
            lp = load_loops_bedpe(f)
        except Exception:
            continue
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

    # (d-f) APA, aggregated per cell type
    by_ct: dict[str, list[np.ndarray]] = {}
    for sid, r in lib.iterrows():
        f = results / f"qc/apa/{sid}.apa.npy"
        if not f.exists():
            continue
        m = np.load(f)
        if m.any():
            by_ct.setdefault(r["cell_type"], []).append(m)

    for i, ct in enumerate(["Naive", "Th17", "Treg"]):
        ax = fig.add_subplot(gs[1, i])
        if i == 0:
            _panel_label(ax, "d")
        mats = by_ct.get(ct)
        if not mats:
            _empty(ax, f"{ct}: no APA")
            continue
        m = np.mean(mats, axis=0)
        win = (m.shape[0] - 1) // 2
        im = ax.imshow(np.log2(m + 1), cmap="Reds", origin="lower",
                       extent=[-win, win, -win, win])
        ax.set_title(f"{ct} (n={len(mats)})", loc="left", color=_colour(ct))
        ax.set_xlabel("bins (10 kb)")
        if i == 0:
            ax.set_ylabel("bins (10 kb)")
        fig.colorbar(im, ax=ax, fraction=0.046, label="log2(1+contacts)")

    fig.suptitle("Figure 3 · Loop calling and aggregate peak analysis", x=0.02,
                 ha="left", fontsize=9, fontweight="bold")
    _save(fig, out)


# ---------------------------------------------------------------- figure 4
def figure4_differential(comparisons: list[str], results: Path, out: Path,
                         fdr: float, lfc: float) -> None:
    n = max(len(comparisons), 1)
    fig = plt.figure(figsize=(2.5 * n, 5.0))
    gs = fig.add_gridspec(2, n, hspace=0.5, wspace=0.4)

    for i, comp in enumerate(comparisons):
        f = results / f"diff/{comp}/differential_loops.tsv"
        ax = fig.add_subplot(gs[0, i])
        if i == 0:
            _panel_label(ax, "a")
        try:
            res = pd.read_csv(f, sep="\t")
        except Exception:
            res = pd.DataFrame()
        if res.empty or "padj" not in res.columns:
            _empty(ax, f"{comp}\nno differential result")
            ax2 = fig.add_subplot(gs[1, i]); _empty(ax2, "")
            continue

        res = res.dropna(subset=["padj", "log2FoldChange"])
        sig = (res["padj"] < fdr) & (res["log2FoldChange"].abs() >= lfc)
        y = -np.log10(res["padj"].clip(lower=np.nextafter(0, 1)))
        ax.scatter(res.loc[~sig, "log2FoldChange"], y[~sig], s=2, c="#d9d9d9",
                   rasterized=True, linewidths=0)
        up = sig & (res["log2FoldChange"] > 0)
        dn = sig & (res["log2FoldChange"] < 0)
        ax.scatter(res.loc[up, "log2FoldChange"], y[up], s=3, c="#C44E52",
                   rasterized=True, linewidths=0, label=f"up ({int(up.sum())})")
        ax.scatter(res.loc[dn, "log2FoldChange"], y[dn], s=3, c="#4C72B0",
                   rasterized=True, linewidths=0, label=f"down ({int(dn.sum())})")
        ax.axhline(-np.log10(fdr), ls="--", lw=0.5, c="black")
        ax.axvline(lfc, ls="--", lw=0.5, c="black")
        ax.axvline(-lfc, ls="--", lw=0.5, c="black")
        ax.set_title(comp.replace("_", " "), loc="left")
        ax.set_xlabel("log2 fold change")
        if i == 0:
            ax.set_ylabel("-log10 adjusted p")
        ax.legend(frameon=False, loc="upper left")

        # span of changed loops — do gained loops act at a different range?
        ax2 = fig.add_subplot(gs[1, i])
        if i == 0:
            _panel_label(ax2, "b")
        if {"start1", "start2"}.issubset(res.columns) and sig.any():
            span = (res["start2"].astype(float) - res["start1"].astype(float)).abs()
            for mask, col, lab in ((up, "#C44E52", "up"), (dn, "#4C72B0", "down")):
                s = span[mask]
                s = s[s > 0]
                if len(s) > 1:
                    ax2.hist(np.log10(s), bins=25, histtype="step", lw=0.9,
                             color=col, density=True, label=lab)
            ax2.set_xlabel("log10 loop span (bp)")
            if i == 0:
                ax2.set_ylabel("density")
            ax2.legend(frameon=False)
        else:
            _empty(ax2, "no significant loops")

    if not comparisons:
        _empty(fig.add_subplot(gs[0, 0]), "no comparisons configured")

    fig.suptitle("Figure 4 · Differential loops", x=0.02, ha="left",
                 fontsize=9, fontweight="bold")
    _save(fig, out)


# ---------------------------------------------------------------- figure 5
def figure5_stripes(lib: pd.DataFrame, results: Path, out: Path) -> None:
    """Stripes, read against the anchor they were called on.

    On CTCF the stripe is directly a loop-extrusion anchor. On H3K27ac the anchors
    are enhancers, extrusion is not what defines them, and fewer/weaker stripes are
    the expected result -- not a failure. The figure is split by mark so the two are
    never averaged together.
    """
    fig = plt.figure(figsize=(7.6, 3.4))
    gs = fig.add_gridspec(1, 3, wspace=0.45)

    stripes: dict[str, pd.DataFrame] = {}
    for sid in lib.index:
        f = results / f"stripes/{sid}/result_filtered.tsv"
        if not f.exists():
            continue
        try:
            df = pd.read_csv(f, sep="\t")
        except Exception:
            continue
        stripes[sid] = df

    # (a) count by sample, coloured by mark
    ax = fig.add_subplot(gs[0, 0]); _panel_label(ax, "a")
    if stripes:
        counts = pd.Series({s: len(d) for s, d in stripes.items()})
        counts = counts.reindex(lib.index).fillna(0).sort_values(ascending=False)
        cols = ["#8172B2" if lib.loc[s, "mark"] == "CTCF" else "#CCB974" for s in counts.index]
        ax.bar(range(len(counts)), counts.values, color=cols,
               edgecolor="black", linewidth=0.3)
        _label_axis(ax, counts.index, axis="x")
        ax.set_ylabel("stripes")
        ax.set_title("Stripe count", loc="left")
        ax.legend(handles=[
            Line2D([0], [0], color="#8172B2", lw=4, label="CTCF anchors"),
            Line2D([0], [0], color="#CCB974", lw=4, label="H3K27ac anchors"),
        ], frameon=False)
    else:
        _empty(ax, "no stripe tables")

    # (b) stripes per mark, normalised for library count
    ax = fig.add_subplot(gs[0, 1]); _panel_label(ax, "b")
    if stripes:
        rec = [{"mark": lib.loc[s, "mark"], "n": len(d)} for s, d in stripes.items() if s in lib.index]
        df = pd.DataFrame(rec)
        if not df.empty:
            for i, (mark, sub) in enumerate(df.groupby("mark")):
                col = "#8172B2" if mark == "CTCF" else "#CCB974"
                ax.scatter(np.full(len(sub), i) + np.random.default_rng(0).normal(0, 0.05, len(sub)),
                           sub["n"], s=16, color=col, edgecolor="black", linewidth=0.3)
                ax.hlines(sub["n"].median(), i - 0.2, i + 0.2, color="black", lw=1)
            ax.set_xticks(range(df["mark"].nunique()))
            ax.set_xticklabels(sorted(df["mark"].unique()))
            ax.set_ylabel("stripes per library")
            ax.set_title("Stripes by anchor type", loc="left")
        else:
            _empty(ax, "no stripe tables")
    else:
        _empty(ax, "no stripe tables")

    # (c) stripe length
    ax = fig.add_subplot(gs[0, 2]); _panel_label(ax, "c")
    drawn = 0
    for sid, d in stripes.items():
        if "length" not in d.columns or d.empty:
            continue
        L = pd.to_numeric(d["length"], errors="coerce").dropna()
        L = L[L > 0]
        if len(L) < 2:
            continue
        col = "#8172B2" if lib.loc[sid, "mark"] == "CTCF" else "#CCB974"
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
    _save(fig, out)


# ---------------------------------------------------------------- main
def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])

    global Q_LABEL
    Q_LABEL = snakemake.params.q_label

    results = Path(snakemake.params.results)
    outdir = Path(snakemake.params.outdir)
    samples = pd.read_csv(snakemake.params.samples_tsv, sep="\t", comment="#").set_index("sample_id")
    min_contacts = int(snakemake.config["hicrep"]["min_contacts_for_scc"])
    fdr = float(snakemake.config["differential"]["fdr"])
    lfc = float(snakemake.config["differential"]["log2fc_min"])
    comparisons = [c["name"] for c in snakemake.config.get("differential", {}).get("comparisons", [])]

    lib = _library_table(samples, results, min_contacts)
    lib.to_csv(snakemake.output.table, sep="\t")
    log.info("library table:\n%s", lib.to_string())

    figure1_library_qc(lib, results, outdir / "figure1_library_qc", min_contacts)
    figure2_reproducibility(lib, results, outdir / "figure2_reproducibility")
    figure3_loops(lib, results, outdir / "figure3_loops_apa")
    figure4_differential(comparisons, results, outdir / "figure4_differential", fdr, lfc)
    figure5_stripes(lib, results, outdir / "figure5_stripes")


# Guarded so the module can be imported by the tests. Snakemake injects
# `snakemake` into the script's globals before executing it.
if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
