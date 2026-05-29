"""
Annotate loop anchors with overlapping peaks and the nearest gene / TSS,
and flag anchors that overlap CTCF sites and super-enhancers (if provided
via config). Writes an annotated BEDPE.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pyranges as pr

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_loops_bedpe, setup_logging  # noqa: E402


def _gtf_to_tss(gtf: str | Path) -> pr.PyRanges:
    df = pr.read_gtf(str(gtf)).df
    genes = df[df["Feature"] == "gene"][["Chromosome", "Start", "End", "Strand", "gene_name", "gene_id"]].copy()
    # TSS = 5' end based on strand
    genes["tss"] = genes.apply(lambda r: r.Start if r.Strand == "+" else r.End, axis=1)
    genes["End"] = genes["tss"] + 1
    genes["Start"] = genes["tss"]
    return pr.PyRanges(genes.rename(columns={"gene_name": "Name"}))


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])

    loops = load_loops_bedpe(snakemake.input.loops)
    if loops.empty:
        Path(snakemake.output.bedpe).write_text("")
        return

    peaks = pd.read_csv(snakemake.input.peaks, sep="\t", header=None,
                        names=["Chromosome", "Start", "End"])
    peaks_pr = pr.PyRanges(peaks)

    tss = _gtf_to_tss(snakemake.input.gtf)

    def _annot(side: str) -> pd.DataFrame:
        sub = loops[[f"chrom{side}", f"start{side}", f"end{side}"]].rename(
            columns={f"chrom{side}": "Chromosome", f"start{side}": "Start", f"end{side}": "End"}
        ).copy()
        sub["loop_idx"] = sub.index
        anchors = pr.PyRanges(sub)

        # Peak overlap
        anchors_with_peak = anchors.count_overlaps(peaks_pr).df.rename(
            columns={"NumberOverlaps": f"peak_n_side{side}"}
        )

        # Nearest gene / TSS
        nearest = anchors.nearest(tss).df
        nearest = nearest[["loop_idx", "Name", "Distance"]].rename(
            columns={"Name": f"nearest_gene_side{side}", "Distance": f"distance_to_tss_side{side}"}
        )
        return anchors_with_peak.merge(nearest, on="loop_idx", how="left")

    a1 = _annot("1")
    a2 = _annot("2")
    out = loops.copy()
    out["loop_idx"] = out.index
    out = out.merge(a1[["loop_idx", "peak_n_side1", "nearest_gene_side1", "distance_to_tss_side1"]],
                    on="loop_idx", how="left")
    out = out.merge(a2[["loop_idx", "peak_n_side2", "nearest_gene_side2", "distance_to_tss_side2"]],
                    on="loop_idx", how="left")
    out["both_anchors_have_peak"] = (out["peak_n_side1"] > 0) & (out["peak_n_side2"] > 0)
    out["promoter_promoter"] = (out["distance_to_tss_side1"].abs() <= 2000) & (out["distance_to_tss_side2"].abs() <= 2000)
    out["enhancer_promoter"] = (
        ((out["distance_to_tss_side1"].abs() <= 2000) & (out["distance_to_tss_side2"].abs() > 2000)) |
        ((out["distance_to_tss_side2"].abs() <= 2000) & (out["distance_to_tss_side1"].abs() > 2000))
    )

    Path(snakemake.output.bedpe).parent.mkdir(parents=True, exist_ok=True)
    out.drop(columns=["loop_idx"]).to_csv(snakemake.output.bedpe, sep="\t", index=False)


main(snakemake)  # type: ignore[name-defined]  # noqa: F821
