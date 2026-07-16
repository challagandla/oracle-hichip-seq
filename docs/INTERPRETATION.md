# Interpreting HiChIP structural tracks

HiChIP is an enrichment assay: its contact map is deliberately concentrated near
the immunoprecipitated histone mark or protein. That makes it excellent for
mark-anchored loops, but it is not a lower-depth substitute for unbiased Hi-C or
Micro-C. Use the following labels and limits when reporting the map-derived tracks.

## Balance convergence comes first

ICE normalization is tested separately at every matrix resolution. Read
`results/qc/balance/<sample>.balance.json` or its TSV companion before comparing a
matrix-derived track:

- `PASS` means cooler recorded `converged=true`; an ICE-balanced value is available.
- `WARN` means the balance attempt did not converge. The raw contact counts remain
  usable, but the attempted weight is removed and must not be interpreted.
- `NOT_ASSESSED` means the convergence evidence is missing. This is unknown, not a
  pass or an observed biological failure.

Expected-cis, insulation, locus heatmaps, virtual 4C, and APA can use a raw-count
fallback, and their tables or plots identify it. That fallback is a coverage-sensitive
description, not an ICE-equivalent normalization. APA under either recorded
normalization is a `DESCRIPTIVE` matched effect with uncertainty, never a universal
pass/fail threshold. E1 has no valid
raw mode in cooltools, so a non-passing 100-kb balance makes E1 `NOT_ASSESSED` rather
than zero. Mustache 1.3.3 likewise requires a converged weight and becomes
`NOT_ASSESSED` at a non-passing configured resolution; it has no raw fallback. These
balance states concern matrix normalization; they do not invalidate
raw-count FitHiChIP, differential, HiCRep, or stripenn results.

## Reported loops and differential hypotheses are different products

The q-filtered, merged FitHiChIP calls are the per-sample loops used for annotation,
APA, locus display, and loop-yield summaries. Differential hypotheses instead come
from FitHiChIP's unthresholded, unmerged all-interaction table. The workflow uses
exact native pixels, primary cis chromosomes, the configured distance range, and a
condition-blind raw-count abundance filter. It does not choose hypotheses using a
q-value, a merged call set, a condition label, or neighbouring-bin tolerance.

For each comparison, read `hypothesis_universe.json` and `design.json` before the
Wald table. The former identifies the source files and independent-filtering rule;
the latter reports the complete-pair count and `analysis_status`. The bundled B2/B3
contrasts have only two complete pairs and are `PILOT_UNDERPOWERED`. Thresholded
rows are exploratory leads, not ordinary publication-strength discoveries. Review
the donor-level `paired_effects.tsv` and Wald uncertainty, and describe any result
as differential mark-associated contact signal because occupancy and 3D contact
frequency are not separated by HiChIP alone.

## Insulation

The 25 kb output is a **local insulation profile**. It can help orient loops and
peaks within a locus and can be compared between biological replicates with the
same mark and adequate usable depth. Do not describe its minima as definitive TAD
boundaries from HiChIP alone. Validate a boundary claim with Hi-C/Micro-C or an
orthogonal boundary assay. When its `normalization` column says `raw-count fallback`,
interpret the profile even more cautiously because coverage and mappability have not
been corrected by converged ICE weights.

## GC-phased E1

GC phasing fixes the otherwise arbitrary sign of E1; it does not remove HiChIP's
anchor enrichment or sequencing-depth effects. Treat the 100 kb output as an
**exploratory A/B-like E1 signal**, compare only same-mark libraries, and confirm
canonical compartment switching with an unbiased contact map. Cross-mark E1
correlations mix biology with antibody target and are not replicate QC.
If `results/qc/compartments/<sample>.cis.eigs.status.json` is `NOT_ASSESSED`, there
is no E1 result to interpret: a header-only table or empty browser track preserves
the workflow schema but contains no measured zero-valued compartments.

## CTCF and H3K27ac are different experiments

Loop and stripe yield depends strongly on usable unique cis contacts. It also
depends on the anchor class: CTCF enriches architectural anchors, whereas H3K27ac
enriches active regulatory elements. A CTCF-versus-H3K27ac difference is therefore
not a biological contrast unless the experiment was explicitly designed with
matched conditions, depth, complexity and replicate structure. The bundled cohort
contains CTCF only in naive cells, so its CTCF stripe output is a technical
demonstration and must not be used as a cross-mark comparison. The libraries remain
visible in QC and raw tables but are excluded from the configured headline stripe
summary.

## Locus plots and virtual 4C

Each configured locus has a named, 0-based viewpoint rather than an inferred window
midpoint. The bundled viewpoints are GENCODE v46 transcription-start coordinates,
and the composite plot includes the same configured GTF as a gene-model track.
Before interpreting a virtual-4C profile, verify adequate support at the viewpoint;
use it to show locus structure, not as a between-sample test without matched depth
and an explicit normalization model. Read the normalization in the plot title and
y-axis. A raw-count fallback is useful for within-map orientation but must not be
presented as ICE-balanced or compared across unequal-depth libraries as if normalized.
