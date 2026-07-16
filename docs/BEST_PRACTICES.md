# Practical review checklist for HiChIP analysis

This checklist describes the safeguards implemented by this repository and the
decisions that still require scientific judgment. The numeric values in
`config/config.yaml` are project screening defaults, not universal pass/fail laws.
Review the underlying counts, assay target, depth, protocol, and replicate design
before accepting or rejecting a library.

For a complete first run, follow [the tutorial](../TUTORIAL.md). For limits on E1,
insulation, stripes, and virtual 4C, read
[the interpretation guide](INTERPRETATION.md).

## Experimental design

- Use independent biological subjects as replicates. Merge technical sequencing
  runs from the same library before analysis; do not count them as independent
  observations.
- Prefer at least three biological replicates per condition for differential work.
  More may be needed for heterogeneous primary tissue or small expected effects.
  Two replicates provide only limited variance estimation and should be treated as a
  pilot.
- Balance library preparation, antibody lot, sequencing, and processing batch across
  conditions. A batch completely confounded with condition cannot be repaired by a
  statistical covariate.
- Keep mark, tissue or cell population, protocol, enzyme, and assembly consistent
  within a differential comparison.
- Add orthogonal assays when the claim requires them. ChIP-seq or CUT&Tag helps
  separate occupancy from contact change; RNA-seq supports expression claims; and
  Hi-C or Micro-C is needed for strong compartment or domain claims.

## Read and contact processing

- Use paired-end data and preserve mate pairing. Appropriate read length and depth
  depend on library insert distribution, genome, assay target, and study aim; inspect
  FastQC, fastp, mapping, and complexity rather than applying one read-length rule.
- Check both mates for gzip integrity, complete FASTQ records, and equal record
  counts. The SRA acquisition rule performs these checks before publishing files.
- Trim adapters and low-quality sequence, then review the amount and reason for
  removal. The configured fastp settings are Q20 and a minimum retained length of
  36 bases; no fixed retention percentage is a universal HiChIP cutoff.
- Align with the configured `bwa-mem2 -SP5M` and parse contacts with pairtools. This
  workflow uses MAPQ 30 and retains `UU` pairs for its high-confidence contact set.
- Assign restriction fragments during contact parsing. The reported dangling-end-
  like, self-circle-like, same-strand-neighbour, and unassigned fractions are
  calculated on deduplicated pairs before the restriction-artifact contact filter;
  they are not a raw-library digestion-efficiency estimate. The analysis contact map
  excludes these neighbouring-fragment artifacts; change that filter only with a
  protocol-specific rationale.
- Deduplicate using both contact ends. Read-level duplicate marking is not a
  substitute for pair-level Hi-C/HiChIP deduplication.

Report QC with explicit denominators:

- raw valid-pair yield = selected deduplicated `UU` contacts / all sequenced read
  pairs from fastp's before-filtering population;
- post-trim valid-pair yield = the same contacts / fastp-retained read pairs
  (descriptive, not the primary gate);
- duplicate rate = duplicate pairs / mapped pairs, using pairtools statistics;
- cis fraction = cis contacts / selected deduplicated `UU` contacts.

The shipped 25% valid-yield, 50% duplicate, and 0.70 cis-fraction gates are cohort
screens. Investigate failures alongside absolute usable cis contacts, complexity,
fragment-orientation QC, and matched libraries before deciding whether to exclude a
sample.

## Contact maps and one-dimensional anchors

- Keep raw counts and balanced values conceptually separate. FitHiChIP and
  differential counting use raw counts; balanced matrices are used for display and
  selected matrix-based QC.
- Remove every contact touching an assembly-blacklist bin before building the
  analysis cooler. Final-loop filtering alone is insufficient because artefacts can
  otherwise influence balancing, coverage-background models, SCC, and stripes.
- Audit ICE convergence independently at every configured resolution. A stored
  `weight` column is not proof of convergence; require cooler's `converged=true`
  metadata before using it.
- Preserve `PASS`, `WARN`, and `NOT_ASSESSED` as distinct balance states. `WARN`
  means the balancing attempt did not converge, while `NOT_ASSESSED` means evidence
  is unavailable. Neither should be silently converted to pass or fail.
- Remove non-passing attempted weights from the published cooler and retain their
  parameters and convergence evidence in the balance JSON/TSV. This prevents tools
  that do not inspect HDF5 attributes from applying a bad weight automatically.
- When balance is not `PASS`, label every permitted raw-count fallback. Expected-cis,
  insulation, locus heatmaps, and virtual 4C remain descriptive but become more
  coverage-sensitive. E1 must instead be `NOT_ASSESSED` because cooltools has no
  valid raw-count `eigs-cis` mode.
- Do not force Mustache through a nonconverged matrix. Mustache 1.3.3 implicitly
  requires `bins/weight`; emit its stable empty TSV and `NOT_ASSESSED` status instead.
- Build all resolutions required by downstream rules from one finest-resolution
  cooler. Do not compare analyses performed on different bin grids as if they were
  the same test.
- Inspect P(s) over the reported distance range and compare like-for-like libraries.
  There is no single acceptable slope for every HiChIP target, cell type, depth, and
  protocol.
- Call one-dimensional anchors from individual read ends after pair-level
  deduplication. The workflow uses MACS3 with `--keep-dup all` because duplicate
  molecules have already been removed at the pair level.
- Match narrow or broad MACS3 mode to the assay target and document changes from the
  configured mark table. Restrict peak calling to the configured primary-assembly
  chromosome view and remove blacklist overlaps before using peaks as anchors.
- The FitHiChIP anchor universe retains bases supported by at least two libraries
  in the configured mark/tissue/protocol/enzyme assay stratum when possible. It is
  pooled across contrast conditions but never across unrelated assay strata. This
  keeps a common search space; it is not proof of within-condition replication.
- Report both read-end FRiP in each sample's own filtered peaks and FRiP in the shared
  assay-stratum consensus anchors. Use the same declared primary-chromosome read
  population for both numerators and the denominator. Treat both as descriptive:
  their magnitudes depend on target, depth, peak width, and how the anchor universe
  was constructed.

## Loop calling and loop QC

- Use the configured FitHiChIP Peak-to-ALL coverage-bias model as the primary caller.
  This workflow calls at 5 kb, 20 kb to 3 Mb, and q < 0.01, retains calls with at
  least six supporting contacts, then filters to cis primary chromosomes and removes
  blacklist overlaps. These are analysis choices, not universal values for every
  HiChIP experiment.
- Keep FitHiChIP's products distinct. The unthresholded, unmerged
  `interactions_FitHiC.bed` table is the source for a differential hypothesis
  universe; the q-filtered `MergeNearContacts` set is for reported calls, APA,
  annotation, and display. The latter can retain multiple representative rows per
  connected neighbourhood, so do not describe it as exactly one row per biological
  loop. Record the exact bin size, distance range, model, FDR threshold, and product
  type with every exported table.
- Use Mustache at 10 kb as independent supporting evidence. Do not silently merge
  callers or imply that caller agreement is a second biological replicate.
- Interpret loop count only with caller-range primary cis depth, anchor opportunity,
  target, distance distribution, and replicate evidence. It is descriptive and is
  not a cohort-independent sample gate.
- Inspect contact-map-held-out APA rather than scoring a map only on its own calls.
  Candidates come from sibling donors and sibling-only anchors. Controls preserve
  chromosome/distance and match anchor class, caller-range marginal visibility,
  blacklist status, and usable matrix coverage.
- Reconcile q-filtered sibling APA calls within the configured one-bin
  reciprocal-anchor tolerance. Keep this separate from differential testing:
  APA tolerates merged-call representative jitter, whereas differential hypotheses
  remain exact native pixels with zero neighbouring-bin tolerance.
- Treat APA as a descriptive matched effect with a loop-bootstrap interval under
  either recorded normalization. It is never a hard library pass/fail. Too little
  matched evidence is `NOT_ASSESSED`, not numeric zero.
- Record that sibling FitHiChIP calls were searched in a shared assay-stratum anchor
  universe that includes the scored sample's peaks. Sibling-only filtering removes
  target-only anchor leakage but does not create strict statistical independence.

## Replicate concordance

- Compare independent donors only within the same configured biological condition,
  tissue, mark, and protocol group. Cross-condition or cross-mark SCC is not
  replicate QC.
- Interpret HiCRep with contact depth. The workflow deterministically downsamples the
  deeper map and marks pairs below `hicrep.min_contacts_for_scc` as depth-confounded.
  That depth is counted on the exact off-diagonal autosomal <=max-distance pixel
  population entering SCC, not on all stored cis contacts.
- Do not use one absolute SCC value as evidence that biological replicates agree.
  Review per-pair SCC, group summaries, depth, P(s), and locus-level concordance
  together.
- The categorical gate uses all depth-qualified pairs: all above threshold is
  `PASS`, all below is `FAIL`, mixed evidence is `DISCORDANT`, and no qualified
  pair is `NOT_ASSESSED`. Minimum, mean, maximum, and group-median SCC are
  descriptive only; selecting the best pair can hide a discordant replicate.
- A low SCC after both libraries clear the depth floor warrants investigation; a low
  SCC driven by a shallow member is not enough to distinguish biological discordance
  from inadequate sequencing.

## Differential contacts

- Define comparisons explicitly. For matched samples, use the pairing variable; the
  bundled cohort fits `~ subject_id + condition` and requires one case and one control
  library per subject.
- Build the testing universe from each selected sample's unthresholded, unmerged
  FitHiChIP all-interaction table. Use exact native 5-kb pixel keys
  (`candidate_tolerance_bins: 0`), primary cis chromosomes, the configured caller
  distance range, and blacklist exclusion. A q-value, merged significant call,
  condition label, or neighbouring-bin tolerance must not decide which hypotheses
  are tested.
- Use complete, explicit paired-subject subsets. The bundled primary comparisons
  select donor2/donor3 (published B2/B3); adding donor1 is a sensitivity analysis
  and must add both conditions for that subject.
- The condition-blind abundance filter requires at least `min_count` raw contacts in
  at least `min_samples` selected libraries. This is an independent-filtering rule,
  not condition-specific significance or replication. Verify
  `hypothesis_universe.json`, `candidate_support.tsv`, and the per-sample normalized
  all-interaction audits.
- Recount every retained pixel from the same blacklist-filtered unbalanced cooler,
  require exact agreement with its source all-interaction count, and require the
  count matrix to contain each expected sample/pixel key exactly once. Missing rows
  must fail rather than silently become zero.
- Filter very low-count candidates before model fitting and verify the exact samples,
  factors, pairing, prefilter, universe contract, and `analysis_status` in
  `design.json`.
- Control the multiple-testing rate and report effect size with uncertainty. The
  shipped reporting threshold is adjusted p < 0.05 and absolute log2 fold change at
  least 1, but the full result table should remain available.
- Treat fewer than three complete pairs as `PILOT_UNDERPOWERED`, not ordinary
  inference. Set `differential.require_publication_ready: true` when a run must fail
  below the configured minimum number of complete pairs. Show per-subject paired
  effects and Wald uncertainty so a headline count cannot hide inconsistent donor
  directions.
- Describe results as differential mark-associated contact signal. HiChIP alone
  cannot determine whether a change came from occupancy, three-dimensional contact
  frequency, or both.
- Do not present a tiny label-permutation space as empirical validation. With three
  matched pairs there are only eight paired assignments; additional donors and
  orthogonal evidence are more informative.

## Figures and reporting

- Show raw valid-ligation yield and usable contact depth beside loop and stripe
  yield, using the configured QC/reporting thresholds rather than hard-coded plot
  cutoffs. Cross-mark or cross-depth bar charts without this context are not
  biological comparisons.
- Include coordinates, resolution, gene models, peaks, and loop arcs in locus plots.
  Use an explicit biological viewpoint for virtual 4C rather than the midpoint of a
  display interval.
- Put the matrix normalization and balance status on every contact heatmap, APA, and
  virtual-4C figure. A reader should not need to inspect a log to learn whether a
  panel is ICE-balanced or a raw-count fallback.
- Treat E1 and insulation from enriched HiChIP maps as exploratory structural context.
  Do not label them as definitive compartments or TAD boundaries without unbiased
  contact data.
- Preserve vector PDF for publication and high-resolution PNG for review. Check axis
  labels, units, colour scales, legends, empty-data panels, and rasterization before
  release.
- Review MultiQC and the per-sample JSON/Markdown reports, including the differential
  design/status section; an aggregate status never replaces inspection of the
  component metrics, denominators, or hypothesis source.
- Confirm that MultiQC exposes the overall balance state and the converged,
  nonconverged, and missing resolutions. Follow any `WARN` or `NOT_ASSESSED` back to
  `results/qc/balance/<sample>.balance.json` or its TSV companion.

## Reproducibility and release hygiene

- Install and verify `environment.runner.yml` plus `workflow/envs/*.yaml` with
  `bash setup.sh --check`. Run the workflow through `bash run.sh` so rules use their
  declared environments.
- Keep deterministic seeds for downsampling, control generation, and bootstrapping.
- Record `config/samples.tsv`, `config/config.yaml`, `config/genome.yaml`, the actual
  reference-file hashes, resolved rule-environment package records, the Git revision,
  and the final provenance manifest with the study.
- Retain the balance JSON/TSV and E1 status JSON with the matrices and reports. They
  are part of the evidence needed to reproduce downstream normalization choices.
- Do not commit FASTQ, BAM, pairs, matrices, results, caches, local notes, credentials,
  or worktrees. Before publishing, inspect `git status` and the staged file list.

## ORACLE prototype export

- Treat the current PyTorch Geometric/HDF5 output as a prototype structural export,
  not a complete multimodal model input. Its node channels are peak-overlap density,
  insulation, and E1; continuous measured mark signal is not included.
- Keep the `.pt`, `.h5`, and per-sample manifest together. The manifest hashes the
  declared export inputs and records feature/edge channel order and limitations.
- Read `node_feature_availability`, per-node `x_observed_mask`, and `blacklist_mask`
  before using fixed channels. Insulation records ICE versus raw fallback; E1 is
  unavailable when its required 100-kb balance did not pass.
- Use the separate PyG/HDF5 `loop` and `adjacent` relations. Coarse duplicate loop
  pairs are consolidated with `fine_loop_count`, and `contained_by`/`contains`
  relations connect successive resolutions.
- Do not infer enhancer targets from nearest-gene annotation alone, and do not treat
  zero-filled missing structural features as measured zero signal.
- Read [the ORACLE export contract](ORACLE_INTEGRATION.md) before downstream use.

## Method references

- [HiChIP assay](https://doi.org/10.1038/nmeth.3999)
- [FitHiChIP method](https://doi.org/10.1038/s41467-019-11950-y) and
  [documentation](https://ay-lab.github.io/FitHiChIP/html/index.html)
- [pairtools documentation](https://pairtools.readthedocs.io/en/latest/)
- [cooler](https://open2c.github.io/cooler/) and
  [cooltools](https://cooltools.readthedocs.io/en/latest/)
- [ENCODE histone ChIP-seq standards](https://www.encodeproject.org/chip-seq/histone/)
