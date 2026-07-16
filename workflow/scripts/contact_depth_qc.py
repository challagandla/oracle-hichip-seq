"""Report contact-depth populations that match downstream caller search spaces."""
import sys
from pathlib import Path

import cooler
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from utils import setup_logging, write_json  # noqa: E402


def contact_populations(
    clr: cooler.Cooler,
    chromosomes: list[str],
    lower_distance: int,
    upper_distance: int,
) -> dict[str, int]:
    """Count stored contacts without mixing trans/diagonal/out-of-range pixels."""
    if lower_distance < 0 or upper_distance < lower_distance:
        raise ValueError("contact-distance bounds must satisfy 0 <= lower <= upper")
    selector = clr.matrix(balance=False, as_pixels=True, join=False)
    primary_cis_offdiagonal = 0
    caller_space = 0
    for chrom in chromosomes:
        table = selector.fetch(chrom)
        offsets_bp = (
            table["bin2_id"].to_numpy(dtype=np.int64)
            - table["bin1_id"].to_numpy(dtype=np.int64)
        ) * int(clr.binsize)
        counts = table["count"].to_numpy(dtype=np.int64)
        off_diagonal = offsets_bp > 0
        primary_cis_offdiagonal += int(counts[off_diagonal].sum())
        in_caller_space = (
            off_diagonal
            & (offsets_bp >= int(lower_distance))
            & (offsets_bp <= int(upper_distance))
        )
        caller_space += int(counts[in_caller_space].sum())
    return {
        "primary_cis_offdiagonal_contacts": primary_cis_offdiagonal,
        "fithichip_distance_range_contacts": caller_space,
    }


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    resolution = int(snakemake.params.resolution)
    clr = cooler.Cooler(
        f"{snakemake.input.mcool}::resolutions/{resolution}"
    )
    chromosomes = [
        line.split("\t")[0]
        for line in Path(snakemake.input.view).read_text().splitlines()
        if line.strip()
    ]
    counts = contact_populations(
        clr,
        chromosomes,
        int(snakemake.params.lower_distance),
        int(snakemake.params.upper_distance),
    )
    write_json(
        {
            "schema": "oracle-hichip-contact-depth-v1",
            "sample": snakemake.wildcards.sample,
            "resolution_bp": resolution,
            "chromosomes": chromosomes,
            "fithichip_distance_range_bp": [
                int(snakemake.params.lower_distance),
                int(snakemake.params.upper_distance),
            ],
            **counts,
        },
        snakemake.output.json,
    )


if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
