# ORACLE HiChIP-seq

A reproducible, end-to-end HiChIP workflow from paired FASTQ or SRA accessions to
quality-controlled contact maps, assay-stratum consensus anchors, filtered FitHiChIP
loops, donor-aware differential contacts, locus figures, MultiQC, and ORACLE graph
exports.

HiChIP enriches chromatin contacts associated with an immunoprecipitated protein or
histone mark. It is excellent for questions such as “which promoters contact active
H3K27ac enhancers?” but it is not unbiased Hi-C. This workflow therefore treats
loops as the primary result and labels insulation and E1 as exploratory context.

New users should follow the [complete tutorial](TUTORIAL.md). Scientific limits and
safe language are summarized in [docs/INTERPRETATION.md](docs/INTERPRETATION.md).

## What the workflow does

| Stage | Main method | Important safeguard |
|---|---|---|
| Read acquisition/QC | SRA Toolkit, FastQC, fastp | Paired FASTQ integrity and atomic publication |
| Alignment | bwa-mem2 `-SP5M` | MAPQ 30; Hi-C-aware alignment flags |
| Contact parsing | pairtools | Pair-level deduplication, UU selection, fragment-artifact filtering |
| Contact maps | cooler/cooltools | Blacklist-touching contacts removed before analysis; only converged ICE weights are used, with explicit fallback |
| Anchors | MACS3 | Read ends, primary-assembly contigs, `--keep-dup all`, blacklist removal |
| Anchor consensus | configurable assay strata | Shared across contrast conditions, separated by mark/tissue/protocol/enzyme |
| Loops | FitHiChIP 11.0 Peak-to-ALL | Coverage-bias model, checksum-pinned source, primary chromosomes, blacklist filter |
| Independent cross-check | Mustache at 10 kb | Secondary evidence; never replaces FitHiChIP |
| Loop QC | contact-map-held-out APA, HiCRep | Sibling-only anchors; anchor/visibility-matched controls; deterministic downsampling |
| Differential contacts | FitHiChIP all-interaction tables, pyDESeq2 | Exact native pixels, condition-blind abundance filter; `~ subject_id + condition` |
| Structural context | insulation, GC-phased E1, stripenn | Explicit HiChIP interpretation limits |
| Reporting/export | MultiQC, PDF/PNG, PyG/HDF5 | Correct denominators, hierarchical graph schema, run manifest |

The bundled example is the public GSE101498 primary CD4 T-cell cohort: three
H3K27ac biological replicates each for naive, Th17, and Treg cells, plus two naive
CTCF libraries. Technical sequencing runs are merged within a biological library;
donors remain separate. Every library is processed and shown in QC. The primary
differential contrasts explicitly use the complete B2/B3 donor pairs, matching the
published differential-table design; the shallow B1 pair is retained as QC evidence
and can be added only as a declared complete-pair sensitivity analysis. Two pairs
are not enough for publication-strength inference: the bundled contrasts are
labelled `PILOT_UNDERPOWERED`, and their model output is exploratory.

## Requirements

- Linux or WSL2 with Bash
- Conda, Mamba, or Miniforge (the installer can install Miniforge)
- At least 8 CPU cores for a practical run; 24–32 recommended
- About 20 GB for software environments
- Storage appropriate to the cohort. For the complete bundled public cohort,
  reserve roughly 1 TB of working space; retained outputs can occupy hundreds of
  gigabytes. A small pilot needs much less.

No global Python or R packages are required.

## Quick start

```bash
git clone https://github.com/challagandla/oracle-hichip-seq.git
cd oracle-hichip-seq

# Install the runner, 19 rule environments, and the checksum-pinned FitHiChIP source.
bash setup.sh

# Download/index hg38 references. Run this after setup so cooler/bwa are available.
bash prepare_references.sh hg38

# Inspect the complete default DAG before doing work.
bash run.sh --dry-run --cores 1

# Run the bundled public cohort.
bash run.sh --cores 24
```

`bash run.sh` activates nothing in the parent shell. Snakemake selects an isolated
environment for each rule. `rule all` is the explicit default target, so the final
command means the complete workflow, not the first rule in the Snakefile.

Before a large run, read the tutorial and edit:

- `config/samples.tsv` — one row per biological library
- `config/config.yaml` — analysis choices and comparisons
- `config/genome.yaml` — reference paths and assembly

Useful checks:

```bash
bash setup.sh --check
bash run.sh --dry-run --summary
bash run.sh --cores 4 multiqc
bash test.sh
```

## Core analysis contracts

### Valid-pair and duplicate metrics

- Raw valid-pair yield = final deduplicated valid-ligation UU contacts / all
  sequenced read pairs from fastp's before-filtering population.
- Post-trim valid-pair yield uses fastp-retained read pairs and is descriptive.
- Duplicate rate = pairtools duplicate pairs / mapped pairs.
- Cis fraction = cis contacts / selected deduplicated valid-ligation UU contacts.

These populations must not be mixed. The per-sample QC JSON includes the underlying
input, mapped, duplicate, final high-confidence, and cis counts, not only percentages.

### Balance convergence and fallback

ICE balancing is attempted and audited separately at every configured matrix
resolution. A stored `weight` column is not enough: the workflow checks cooler's
convergence metadata and writes both a machine-readable report and a table:

```text
results/qc/balance/<sample>.balance.json
results/qc/balance/<sample>.balance.tsv
```

- `PASS` means that resolution has a weight with `converged=true` and balanced
  values may be used.
- `WARN` means balancing was attempted but did not converge. The raw contact map is
  still valid, but the attempted weights are not published or used.
- `NOT_ASSESSED` means a weight or its convergence evidence is missing. It must not
  be treated as a pass.

The sample-level status summarizes all configured resolutions. It is a matrix
normalization status, not an automatic reason to discard a library. FitHiChIP,
differential counting, HiCRep, and stripenn continue to use their intended raw-count
inputs. Matrix-derived outputs fail safely: expected-cis, insulation, locus heatmaps,
virtual 4C, and APA can use a clearly labelled raw-count fallback. APA is always a
descriptive matched effect size rather than a universal sample gate; normalization
is recorded alongside it.
Mustache 1.3.3 also has no valid raw mode: if its configured 10-kb balance is not
`PASS`, the secondary cross-check is `NOT_ASSESSED` and emits a header-only TSV plus
status JSON. When it does run, a deterministic manager-free threaded adapter avoids
Mustache's otherwise unnecessary local socket server on restricted cluster nodes;
the status records the backend and output validation. Cooltools E1 has no valid raw
mode, so a non-passing 100-kb balance produces an empty
schema-valid E1 table and `NOT_ASSESSED`, never invented zero-valued E1 measurements.
Normalization and balance status are carried into tables, plots, MultiQC, provenance,
and ORACLE feature-availability metadata.

### Anchors and loops

Contacts touching an assembly-blacklist bin are removed before the base cooler is
built, so they cannot influence balancing, FitHiChIP background estimation, HiCRep,
stripes, or differential counts. MACS3 receives individual read ends from the
deduplicated contact library. Pair-level deduplication has already removed duplicate
molecules, so MACS3 is told to keep all remaining tags. Peak calling is restricted
to the configured primary-assembly chromosome view and blacklisted peaks are
removed. FitHiChIP uses one consensus anchor
set per configured assay stratum (`mark`, `tissue`, `library_protocol`, and
`restriction_enzyme` by default), pooled across contrast conditions. Unrelated
tissues or protocols therefore do not share a search space, while case and control
do not define different peak universes.

FitHiChIP receives the unbalanced 5-kb cooler, not a HiC-Pro validPairs matrix. Its
Peak-to-ALL coverage-bias model produces two deliberately separate products:

- an unthresholded, unmerged all-interaction table used to define a common
  differential-testing universe without selecting on q-value or condition; and
- q-filtered, merged calls used for per-sample loop reporting, annotation, APA, and
  visualization.

Both products are normalized to explicit schemas and audited. Reported calls are
limited to primary autosomes and chrX, exclude blacklist overlaps, and meet the
configured minimum contact count. Per-sample anchor FRiP and shared-consensus-anchor
FRiP are both reported against the same primary-chromosome read-end denominator;
neither is a universal pass/fail threshold.

### Replicate and APA QC

HiCRep compares biological donors only within the configured condition, tissue,
mark, and protocol group. Its stochastic depth matching uses a stable pair-specific
seed. Low-depth comparisons are reported as depth-confounded rather than forced to
pass or fail. Among depth-qualified pairs,
a sample and its replicate group pass only when every relevant pair clears the SCC
threshold, fail when every pair is below it, and are `DISCORDANT` when results are
mixed. Minimum, mean, best, and group-median SCC values are descriptive summaries,
not alternative gates.

APA is held out at the scored-contact-map level: candidates come from sibling donors
and must overlap sibling-only anchors. The sibling FitHiChIP calls were originally
searched in the shared assay-stratum anchor universe, so that residual cohort-level
search-space dependence is recorded rather than hidden. Q-filtered sibling reporting
calls are reconciled within one native FitHiChIP bin per anchor
(`apa.candidate_tolerance_bins: 1`) because merged calls can jitter between donors;
the differential hypothesis universe deliberately remains exact-grid
(`differential.candidate_tolerance_bins: 0`). Random controls preserve
chromosome and distance and match sibling-anchor class, caller-range marginal
visibility, blacklist status, and usable matrix coverage. The matched effect and
loop-bootstrap interval are `DESCRIPTIVE`; insufficient evidence is
`NOT_ASSESSED`. APA never hard-passes or hard-fails a library.

### Differential contacts

Differential hypotheses come from every selected sample's unthresholded FitHiChIP
all-interaction table. Rows are restricted to the exact native 5-kb grid, configured
20 kb–3 Mb search range, primary chromosomes, and non-blacklisted anchors. A pixel
enters the condition-blind common universe only when its raw FitHiChIP count meets
the configured abundance threshold in the configured number of selected libraries.
No q-value, merged significant call, condition label, reciprocal-anchor tolerance,
or transitive neighbourhood is used to choose hypotheses. Counts are then read from
the same blacklist-filtered unbalanced cooler in every sample and must exactly match
the source-table count at the source pixel. Strict key and sample validation prevents
missing rows from becoming invented zeros.

The bundled primary contrasts select `include_subjects: [donor2, donor3]` and fit:

```text
~ subject_id + condition
```

Low-count candidates are filtered before pyDESeq2. `design.json`, the hypothesis
manifest, per-subject paired effects, and the full Wald result table make the tested
population and uncertainty explicit. With only two complete pairs, the bundled
analysis is emitted as `PILOT_UNDERPOWERED`; setting
`differential.require_publication_ready: true` requires at least the configured
number of complete pairs and fails rather than relabelling a pilot. Results are
correctly described as differential mark-associated contact signal: a change can
reflect contact frequency, anchor occupancy, or both.

## Main outputs

```text
results/
├── pairs/                   deduplicated, indexed .pairs.gz
├── cool/                    base .cool and multi-resolution .mcool; only converged ICE weights retained
├── peaks/                   per-library and assay-stratum consensus anchors
├── loops/                   filtered FitHiChIP calls and Mustache cross-checks
├── qc/                      pairtools, balance, restriction, FRiP, P(s), APA, HiCRep, E1/status
├── diff/                    hypothesis manifest, strict counts/design, paired effects, results, MA/volcano
├── stripes/                 stripenn calls and cohort summary
├── viz/                     locus composites and explicit-viewpoint virtual 4C
├── figures/                 five cohort figures as vector PDF and 400-dpi PNG
├── multiqc/                 self-contained MultiQC HTML
├── oracle_cos/              hierarchical PyG/HDF5 v2 graph prototypes + manifests
└── provenance/              config/reference/resolved-environment/output hashes and Git revision
```

Generated data, results, local worktrees, graph indexes, notes, and local development
metadata are ignored by Git. Only source, configuration, tests, and user documentation
belong in the repository.

## Cluster execution

A starter SLURM profile is included:

```bash
bash run.sh --profile profiles/slurm --jobs 100
```

Review memory, runtime, partition, account, and scratch settings before use. Keep the
repository, reference files, Conda cache, and results on storage visible to compute
nodes. `setup.sh` prefetches FitHiChIP so loop jobs do not need internet access.

## What not to claim

- HiChIP E1 is exploratory A/B-like signal, not canonical compartment evidence.
- HiChIP insulation is locus context, not definitive TAD discovery.
- CTCF and H3K27ac loop/stripe yields are not directly comparable in this cohort.
- The two unmatched CTCF libraries remain in QC and raw outputs but are configured
  as demonstrations and excluded from the headline stripe summary.
- Nearest-gene annotation is navigation, not proof of enhancer target assignment.
- A raw-count fallback is not equivalent to ICE balancing; preserve its label and
  interpret coverage-sensitive tracks descriptively.
- A within-study differential result requires biological replication and a valid
  design; technical sequencing runs are not replicates.
- `PILOT_UNDERPOWERED` model output is exploratory screening, not ordinary
  publication-strength statistical evidence, even when a Wald adjusted p-value is
  small.

## References

- Mumbach et al. *Nature Methods* (2016), HiChIP: https://doi.org/10.1038/nmeth.3999
- Mumbach et al. *Nature Genetics* (2017), GSE101498: https://doi.org/10.1038/ng.3963
- Bhattacharyya et al. FitHiChIP: https://doi.org/10.1038/s41467-019-11950-y
- Open2C pairtools: https://pairtools.readthedocs.io/
- cooler/cooltools: https://open2c.github.io/cooler/ and https://cooltools.readthedocs.io/

## License and citation

MIT License. See `LICENSE` and `CITATION.cff`.
