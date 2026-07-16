"""Stream restriction-fragment orientation QC from a restricted .pairs file."""
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from utils import open_text_auto, setup_logging, write_json  # noqa: E402


def classify(chrom1: str, chrom2: str, rfrag1: int, rfrag2: int,
             strand1: str, strand2: str) -> str:
    if rfrag1 < 0 or rfrag2 < 0:
        return "unassigned"
    if chrom1 != chrom2 or abs(rfrag1 - rfrag2) > 1:
        return "regular"
    if strand1 == "+" and strand2 == "-":
        return "dangling_end_like"
    if strand1 == "-" and strand2 == "+":
        return "self_circle_like"
    return "same_strand_neighbour"


def summarise_pairs(path: str | Path) -> dict:
    columns = None
    counts: Counter[str] = Counter()
    with open_text_auto(path) as handle:
        for line in handle:
            if line.startswith("#columns:"):
                columns = line.rstrip().split()[1:]
                continue
            if line.startswith("#"):
                continue
            if columns is None:
                raise ValueError("pairs header lacks a #columns declaration")
            fields = line.rstrip().split("\t")
            row = dict(zip(columns, fields))
            required = {"chrom1", "chrom2", "strand1", "strand2", "rfrag1", "rfrag2"}
            if not required.issubset(row):
                raise ValueError(
                    "pairs file lacks restriction columns; run pairtools restrict before dedup"
                )
            counts[classify(
                row["chrom1"], row["chrom2"], int(row["rfrag1"]), int(row["rfrag2"]),
                row["strand1"], row["strand2"],
            )] += 1
    total = sum(counts.values())
    retained = counts.get("regular", 0)
    return {
        "population": "post_dedup_pre_contact_filter_UU_pairs",
        "denominator_description": (
            "Deduplicated high-confidence UU pairs after restriction-fragment "
            "assignment and before restriction-artifact contact filtering"
        ),
        "total_deduplicated_uu_pairs": total,
        "valid_ligation_pairs": retained,
        "restriction_artifact_pairs": total - retained,
        "valid_ligation_fraction": retained / total if total else 0.0,
        "counts": dict(sorted(counts.items())),
        "fractions": {k: v / total if total else 0.0 for k, v in sorted(counts.items())},
        "interpretation": (
            "Regular pairs are retained for contact matrices. Neighbouring-fragment "
            "dangling-end-like, self-circle-like, same-strand, and unassigned pairs "
            "are reported here and excluded by the default contact-map filter."
        ),
    }


def main(snakemake) -> None:  # type: ignore[no-untyped-def]
    setup_logging(snakemake.log[0])
    report = {"sample": snakemake.wildcards.sample, **summarise_pairs(snakemake.input.pairs)}
    write_json(report, snakemake.output.json)


if "snakemake" in globals():
    main(snakemake)  # type: ignore[name-defined]  # noqa: F821
