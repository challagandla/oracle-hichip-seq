"""
Build a pyGenomeTracks composite figure for a region of interest:
  - HiChIP contact matrix heatmap
  - exploratory local insulation track
  - 1D MACS3 peaks
  - arc plot of FitHiChIP loops
  - GENCODE gene models
"""
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from balance_utils import load_balance_report, resolution_balance  # noqa: E402
from utils import load_loops_bedpe, select_insulation_column, setup_logging  # noqa: E402

LOOPS_TEMPLATE = """\
[loops]
file = {loops}
title = FitHiChIP loops
height = 2
file_type = links
links_type = arcs
line_width = 0.5
color = #c0392b

[spacer]
height = 0.05
"""

INI_TEMPLATE = """\
[x-axis]

[spacer]
height = 0.1

[matrix]
file = {mcool}::resolutions/{res}
title = {matrix_title}
height = 5
file_type = hic_matrix
colormap = Reds
transform = log1p
depth = {depth}
show_masked_bins = false

[spacer]
height = 0.05

[insulation]
file = {insul_bdg}
title = {insulation_title}
height = 1.5
color = #444
file_type = bedgraph
type = line

[spacer]
height = 0.05

[peaks]
file = {peaks}
title = peaks ({mark})
height = 0.8
color = #1f5fbf
file_type = bed
display = collapsed

{loops_section}

[genes]
file = {gtf}
title = GENCODE genes
height = 3.5
file_type = gtf
prefered_name = gene_name
merge_transcripts = true
merge_overlapping_exons = true
labels = true
fontsize = 7
style = UCSC
"""


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])

    sample = snakemake.wildcards.sample
    region = snakemake.params.region
    mark = snakemake.params.mark
    res = int(snakemake.params.res)
    region_span = int(region["end"]) - int(region["start"])
    if region_span <= 0:
        raise ValueError(f"invalid visualisation region {region!r}")
    matrix_depth = min(1_500_000, region_span)

    import pandas as pd
    insul = pd.read_csv(snakemake.input.insul, sep="\t")
    col = select_insulation_column(insul)
    balance_report = load_balance_report(snakemake.input.balance)
    matrix_balance = resolution_balance(balance_report, res)
    if "normalization" not in insul.columns or "balance_status" not in insul.columns:
        raise ValueError("insulation table lacks balance normalization annotations")
    insulation_normalizations = set(insul["normalization"].dropna().astype(str))
    insulation_statuses = set(insul["balance_status"].dropna().astype(str))
    if len(insulation_normalizations) != 1 or len(insulation_statuses) != 1:
        raise ValueError("insulation table has inconsistent balance annotations")
    insulation_normalization = next(iter(insulation_normalizations))
    insulation_status = next(iter(insulation_statuses))
    # Region-specific sidecars avoid same-sample plot jobs racing while writing
    # the same temporary BED/BEDPE files in parallel.
    bedgraph = Path(snakemake.output.insulation_bedgraph)
    bedgraph.parent.mkdir(parents=True, exist_ok=True)
    local_insul = insul.loc[
        (insul["chrom"] == region["chrom"])
        & (insul["start"] < int(region["end"]))
        & (insul["end"] > int(region["start"])),
        ["chrom", "start", "end", col],
    ].dropna()
    if local_insul.empty:
        raise ValueError(
            f"no insulation values overlap {region['chrom']}:"
            f"{region['start']}-{region['end']}"
        )
    local_insul.to_csv(bedgraph, sep="\t", header=False, index=False)

    # FitHiChIP's first row is a plain-text header, not a BED comment. Passing the
    # raw file to pyGenomeTracks makes its bedtools prefilter fail and can turn a
    # large loop set into a slow, warning-filled full-file scan. Emit strict BEDPE
    # coordinates for plotting and omit the arc track cleanly when no loops exist.
    loop_df = load_loops_bedpe(snakemake.input.loops)
    if not loop_df.empty:
        same_chrom = (
            (loop_df["chrom1"] == region["chrom"])
            & (loop_df["chrom2"] == region["chrom"])
        )
        left = loop_df[["start1", "start2"]].min(axis=1)
        right = loop_df[["end1", "end2"]].max(axis=1)
        overlaps_view = (left < int(region["end"])) & (right > int(region["start"]))
        # pyGenomeTracks intentionally omits an arc whose two anchors both lie
        # outside the view and whose span merely crosses the entire panel.
        surrounds_view = (left < int(region["start"])) & (right > int(region["end"]))
        loop_df = loop_df.loc[same_chrom & overlaps_view & ~surrounds_view]

    loop_bedpe = Path(snakemake.output.loop_bedpe)
    loop_bedpe.parent.mkdir(parents=True, exist_ok=True)
    loop_df[["chrom1", "start1", "end1", "chrom2", "start2", "end2"]].to_csv(
        loop_bedpe, sep="\t", header=False, index=False
    )
    if loop_df.empty:
        loops_section = ""
    else:
        loops_section = LOOPS_TEMPLATE.format(loops=loop_bedpe)

    ini_text = INI_TEMPLATE.format(
        mcool=snakemake.input.mcool, res=res, depth=matrix_depth, sample=sample,
        matrix_title=(
            f"HiChIP {sample} — {matrix_balance['normalization']} "
            f"(balance {matrix_balance['status']})"
        ),
        insulation_title=(
            "local insulation — "
            f"{insulation_normalization} (balance {insulation_status}; "
            "exploratory in HiChIP)"
        ),
        peaks=snakemake.input.peaks, mark=mark,
        loops_section=loops_section, insul_bdg=str(bedgraph),
        gtf=snakemake.input.gtf,
    )
    Path(snakemake.output.ini).write_text(ini_text)

    coord = f"{region['chrom']}:{region['start']}-{region['end']}"
    cmd = [
        "pyGenomeTracks", "--tracks", str(snakemake.output.ini),
        "--region", coord,
        "--outFileName", str(snakemake.output.png),
        "--dpi", "200", "--width", "20", "--height", "12",
        "--trackLabelFraction", "0.12",
    ]
    subprocess.run(cmd, check=True)


# Guarded so the module can be imported by the tests. Snakemake injects
# `snakemake` into the script's globals before executing it.
if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
