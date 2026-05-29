"""
Build a pyGenomeTracks composite figure for a region of interest:
  - HiChIP contact matrix heatmap
  - insulation score track
  - 1D MACS2 peaks
  - arc plot of FitHiChIP loops
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils import setup_logging  # noqa: E402

INI_TEMPLATE = """\
[x-axis]

[spacer]
height = 0.1

[matrix]
file = {mcool}::resolutions/{res}
title = HiChIP {sample}
height = 5
file_type = hic_matrix
colormap = Reds
transform = log1p
depth = 1500000
show_masked_bins = false

[spacer]
height = 0.05

[insulation]
file = {insul_bw}
title = insulation
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

[loops]
file = {loops}
title = FitHiChIP loops
height = 2
file_type = links
links_type = arcs
line_width = 0.5
color = #c0392b
"""


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])

    sample = snakemake.wildcards.sample
    region = snakemake.params.region
    res = int(snakemake.params.res)
    # quick bedgraph of insulation for pgt
    import pandas as pd
    insul = pd.read_csv(snakemake.input.insul, sep="\t")
    if "log2_insulation_score" in insul.columns:
        col = "log2_insulation_score"
    elif "insulation" in insul.columns:
        col = "insulation"
    else:
        col = insul.columns[-1]
    bedgraph = Path(snakemake.output.ini).parent / f"{sample}.insulation.bdg"
    insul[["chrom", "start", "end", col]].dropna().to_csv(bedgraph, sep="\t", header=False, index=False)

    # find sample's mark from the samples sheet (passed via Snakemake config)
    mark = "peak"  # default label; sample sheet not in scope here

    ini_text = INI_TEMPLATE.format(
        mcool=snakemake.input.mcool, res=res, sample=sample,
        peaks=snakemake.input.peaks, mark=mark,
        loops=snakemake.input.loops, insul_bw=str(bedgraph),
    )
    Path(snakemake.output.ini).write_text(ini_text)

    coord = f"{region['chrom']}:{region['start']}-{region['end']}"
    cmd = [
        "pyGenomeTracks", "--tracks", str(snakemake.output.ini),
        "--region", coord,
        "--outFileName", str(snakemake.output.png),
        "--dpi", "200", "--width", "20", "--height", "12",
    ]
    subprocess.run(cmd, check=True)


main(snakemake)  # type: ignore[name-defined]  # noqa: F821
