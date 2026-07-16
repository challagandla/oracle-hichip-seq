# HiChIP from reads to defensible results

This tutorial is written for a first-time user, but it does not hide the decisions
that matter in a publication-quality analysis. It explains what HiChIP measures,
how to install and run the workflow, how each stage works, how to judge quality,
and which conclusions the data can and cannot support.

The short version is:

```bash
bash setup.sh
bash prepare_references.sh hg38
bash run.sh --dry-run --cores 1
bash run.sh --cores 24
```

Do not start the last command until you have read the experimental-design and
configuration sections below. The bundled cohort is real and large.

## 1. What HiChIP measures

Hi-C measures chromatin contacts genome-wide. ChIP-seq enriches DNA associated
with a protein or histone modification. HiChIP combines these ideas: proximity-
ligated DNA is immunoprecipitated, so the final contact map is enriched for
contacts associated with the chosen mark.

For H3K27ac HiChIP, the natural result is a set of active-regulatory contacts such
as enhancer-promoter loops. For CTCF HiChIP, architectural loops and stripes are
more direct targets. The assay is not an unbiased low-depth Hi-C experiment:
insulation and E1 can provide locus context, but canonical TAD or compartment
claims should be confirmed with Hi-C or Micro-C.

Three quantities are easy to confuse:

- Anchor enrichment: how strongly the mark is present at a locus.
- Contact frequency: how many ligation molecules connect two loci.
- Loop significance: whether a contact is enriched relative to the distance and
  coverage background modeled by the caller.

A differential H3K27ac HiChIP contact can change because anchor occupancy changed,
because 3D contact frequency changed, or because both changed. Report it as
“differential H3K27ac-associated contact signal” unless orthogonal data separate
those mechanisms.

## 2. Design the experiment before running software

### Biological and technical replicates

A biological replicate is an independently collected donor, animal, culture, or
experimental unit. A technical replicate is another sequencing run from the same
library.

- Keep biological replicates as separate rows in `config/samples.tsv`.
- Put comma-separated SRA run accessions from the same library in one row. The
  workflow merges those runs before analysis.
- Never count technical runs as independent n in a differential test.
- Aim for at least three biological replicates per condition. Two complete paired
  units can be screened, but the workflow labels that result
  `PILOT_UNDERPOWERED`; it is not publication-strength inference.

### Match the design

Conditions in a differential comparison should use the same:

- immunoprecipitated mark or protein;
- tissue/cell population definition;
- library protocol and restriction enzyme;
- reference assembly;
- processing batch, or a non-confounded batch design that can be modeled.

The workflow refuses comparisons that mix mark, tissue, or protocol. Multiple
batches are accepted only when `batch` is listed as a non-confounded covariate and
the complete design is full rank. The bundled GSE101498 H3K27ac conditions are
matched by donor, so the model is `~ subject_id + condition`.

This is a reanalysis, not an exact reproduction of every published contrast. The
study reports Naive, Th17, and Treg cells from three healthy subjects, and the
sample numbering supports donor pairing. Its published differential-loop table used
subjects B2/B3. The bundled workflow processes and reports all three donors, but its
primary differential comparisons explicitly select the complete B2/B3 pairs because
B1—especially Naive B1—is much shallower. Adding B1 is an optional sensitivity
analysis: include the complete case/control pair, document the design change, and do
not silently add only the better member of a donor pair. The default B2/B3 contrasts
are intentionally pilot analyses. They preserve the published subset and paired
design, but two complete pairs do not justify ordinary confirmatory language.

### Controls and orthogonal data

Useful additions include an input/IgG control for 1D enrichment, ChIP-seq or CUT&Tag
for occupancy, RNA-seq for expression, and unbiased Hi-C/Micro-C for structural
claims. This workflow can run without those assays, but it will not pretend their
questions have been answered.

## 3. Plan compute and storage

For a small pilot, 8 cores, 32 GB RAM, and tens of gigabytes may be sufficient. A
deep multi-library cohort should use a workstation or cluster.

For the full bundled 11-library public cohort, a practical plan is:

- 24–32 CPU cores;
- 128 GB RAM available for the largest concurrent tasks;
- roughly 1 TB of working storage, preferably fast local or parallel scratch;
- about 20 GB for the runner and 19 isolated rule environments.

Retained outputs can occupy hundreds of gigabytes, while alignment, sorting, SRA
conversion, and temporary files raise the peak requirement further. Do not place
the run on a nearly full home directory.

Check space before starting:

```bash
df -h .
du -sh .snakemake 2>/dev/null || true
```

## 4. Install the workflow

Clone the repository and enter it:

```bash
git clone https://github.com/challagandla/oracle-hichip-seq.git
cd oracle-hichip-seq
```

Install:

```bash
bash setup.sh
```

The installer:

1. finds Conda/Mamba, or installs a checksum-verified Miniforge without modifying
   shell startup files;
2. creates the small `oracle-hichip-runner` environment;
3. creates 19 rule-specific environments;
4. runs a real import/executable smoke check in each environment;
5. downloads and verifies the pinned FitHiChIP 11.0 source archive for offline
   compute-node use.

You do not need to run `conda activate`. `run.sh` selects the runner and Snakemake
selects each rule environment.

Verify later with:

```bash
bash setup.sh --check
```

Useful installer modes:

```bash
bash setup.sh --runner-only   # small runner now; rule envs on first execution
bash setup.sh --help
```

## 5. Prepare references

For the bundled human cohort:

```bash
bash prepare_references.sh hg38
```

For mouse:

```bash
bash prepare_references.sh mm10
```

The helper downloads versioned FASTA/annotation files, verifies official MD5 or
SHA-256 checksums, tests gzip integrity, decompresses atomically, builds samtools
and bwa/bwa-mem2 indexes, writes chromosome sizes, and creates an MboI digest BED.
The blacklist is pinned to an immutable repository commit.

Reference files live under `resources/<assembly>/` and are ignored by Git. Check:

```bash
test -s resources/hg38/GRCh38.primary_assembly.genome.fa
test -s resources/hg38/hg38.chrom.sizes
test -s resources/hg38/MboI.digest.hg38.bed.gz
test -s resources/hg38/gencode.v46.primary_assembly.annotation.gtf.gz
```

T2T-CHM13 is present as a configuration template but is not downloaded by the
helper. Supply its FASTA, annotation, blacklist, chromosome sizes, indexes, and
MboI digest at the paths in `config/genome.yaml` before selecting that assembly.

## 6. Describe samples correctly

Open `config/samples.tsv`. One row is one biological library.

| Column | Meaning |
|---|---|
| `sample_id` | Unique path-safe name using letters, numbers, dot, dash, underscore |
| `srr` | Optional SRA run(s); comma-separate technical runs from the same library |
| `subject_id` | Donor/animal/experimental-unit identifier used for pairing |
| `cell_type` | Biological condition or cell population |
| `tissue` | Tissue class used to prevent invalid comparisons |
| `disease` | Disease/phenotype metadata |
| `replicate` | Biological replicate number |
| `mark` | H3K27ac, CTCF, or another mark configured under `macs3.marks` |
| `fastq_r1`, `fastq_r2` | Paired FASTQ paths, whether local or generated from SRA |
| `batch` | Library/preparation batch |
| `library_protocol` | Protocol and enzyme, for example `HiChIP_MboI` |
| `restriction_enzyme` | Exact enzyme used for the library; must match the configured digest |
| `notes` | Human-readable context; not used as a design factor |

For local reads, leave `srr` empty and use paths such as:

```text
data/fastq/my_sample_R1.fastq.gz
data/fastq/my_sample_R2.fastq.gz
```

For SRA, the final paths must be `<fastq_dir>/<sample_id>_R1.fastq.gz` and
`<fastq_dir>/<sample_id>_R2.fastq.gz`; the preflight validator checks this contract.
Enter the accession(s) in `srr`. The fetch rule downloads into a temporary directory, checks that both mates
exist, verifies gzip and FASTQ record structure, compares R1/R2 counts, and only
then publishes the final files.

## 7. Configure the analysis

The main settings are in `config/config.yaml`.

### Assembly and outputs

```yaml
assembly: hg38
samples_tsv: config/samples.tsv
genome_yaml: config/genome.yaml
results_dir: results
fastq_dir: data/fastq
```

The selected genome entry supplies the MACS3 effective-genome label and restriction
digest. `restriction_enzyme` in every sample row must match it. The bundled
visualization loci are hg38-specific; replace the loci and set `viz.assembly` when
using another assembly.

### Contact parsing

```yaml
pairtools:
  min_mapq: 30
  walks_policy: mask
  keep_pair_types: [UU]
  filter_restriction_artifacts: true
```

`UU` means both ends map uniquely. Pairtools assigns MboI fragments before
deduplication and reports dangling-end-like, self-circle-like,
same-strand-neighbour, regular, and unassigned fractions. By default, neighbouring-
fragment artifacts and unassigned pairs are removed from contact maps and loop
calling. The pre-filter deduplicated UU alignments still feed 1D anchor calling.

The explicit fastp adapter sequences are the bundled TruSeq defaults. Confirm them
against the actual library kit before running a new cohort; change both adapter keys
when a different construct was used.

### Contact-map resolutions

```yaml
cooler:
  bin_sizes_kb: [5, 10, 25, 50, 100, 250, 500, 1000, 2500]
  oracle_bin_sizes_kb: [5, 25, 100, 1000]
```

The preflight validator checks that every resolution required by FitHiChIP, APA,
HiCRep, insulation, E1, Mustache, stripenn, and ORACLE exists.

Balancing is also assessed per resolution, because a high-resolution matrix can be
too sparse for ICE even when a coarser matrix converges. Do not infer success from
the presence of a `weight` column. The workflow checks cooler's convergence
attribute, removes non-passing weights from the published `.mcool`, and records the
decision for every resolution.

### Replicate grouping

```yaml
hicrep:
  group_by: [cell_type, tissue, disease, cancer_type, mark, library_protocol, restriction_enzyme]
```

These columns define which libraries are biological replicates for HiCRep. Keep
`subject_id`, `replicate`, and `batch` out of this list: different donors are the
replicates, and a batch difference should remain visible rather than splitting the
comparison. Add an explicit biological condition column when the sample sheet needs
to distinguish time points, treatments, or phenotypes not captured above.

### Consensus-anchor strata

```yaml
anchor_consensus:
  group_by: [mark, tissue, library_protocol, restriction_enzyme]
```

This is an assay-stratum key, not a biological contrast. Peaks are pooled across
case/control conditions inside the stratum, preventing each arm from defining a
different FitHiChIP hypothesis universe. Unrelated tissues, marks, protocols, or
enzymes remain separate. Preflight rejects a comparison whose changing filter
column is also placed in this key.

### FitHiChIP

```yaml
fithichip:
  bin_size: 5000
  lower_distance: 20000
  upper_distance: 3000000
  fdr_threshold: 0.01
  interaction_type: Peak-to-ALL
  background_type: Coverage_Bias
  merge_nearby: true
  min_reads: 6
```

FitHiChIP receives the blacklist-filtered unbalanced cooler and the consensus peaks
for that sample's configured assay stratum. The final file is FitHiChIP's
`MergeNearContacts` call set, normalized to a standard BEDPE-like header and
filtered to primary chromosomes, non-blacklisted anchors, and the configured
minimum contact count. That call set can retain several representative rows for
one connected neighbourhood; it is not exactly one row per biological loop.

### Differential comparisons

```yaml
differential:
  paired_by: subject_id
  covariates: []
  hypothesis_source: fithichip_all_interactions
  candidate_tolerance_bins: 0
  min_count: 5
  min_samples: 2
  publication_min_complete_pairs: 3
  require_publication_ready: false
  comparisons:
    - name: Th17_vs_Naive_H3K27ac
      mark: H3K27ac
      include_subjects: [donor2, donor3]
      case_filter:    {cell_type: Th17}
      control_filter: {cell_type: Naive}
```

The two filters select samples from the sheet; `include_subjects` then makes any
published-style or sensitivity subset explicit. Pairing is checked before the DAG
is built and again before pyDESeq2. A covariate that is constant across selected
samples is recorded but omitted from the fitted design; a confounded batch should
be fixed experimentally, not “adjusted” with an impossible model.

`hypothesis_source` and `candidate_tolerance_bins` are scientific contracts, not
tuning conveniences. Differential hypotheses must come from FitHiChIP's
unthresholded all-interaction table at the exact native pixel (`0` tolerance). The
condition-blind abundance filter keeps a pixel when it has at least `min_count`
contacts in at least `min_samples` selected libraries. With fewer than
`publication_min_complete_pairs`, the model may run for exploration but its
`analysis_status` is `PILOT_UNDERPOWERED`. Set `require_publication_ready: true` for
a hard preflight/runtime gate; the bundled two-pair defaults deliberately leave it
`false` so the tutorial pilot can complete without being mislabelled.

Differential comparisons in this release are deliberately paired-only. If your
experiment is unmatched, do not set `paired_by: null` and reuse the paired figures or
power labels; define and validate a separate unmatched design before extending the
workflow. Never remove real pairing merely to make invalid metadata pass.

The bundled `reporting.demonstration_samples` list contains the two unmatched CTCF
libraries. They are still processed and shown in QC and raw tables, but they do not
enter the headline stripe summary. Change this role only when your study has a
matched, adequately powered CTCF design.

### Regions and virtual 4C

Each visualization region needs an explicit 0-based viewpoint:

```yaml
viz:
  assembly: hg38
  regions:
    - name: MYC_locus
      chrom: chr8
      start: 127235434
      end: 128242951
      viewpoint: 127735433
      viewpoint_label: MYC TSS
```

The plotting window and biological viewpoint are different concepts. The workflow
will not silently use the midpoint of a one-megabase display interval.

## 8. Validate before computing

Run a dry run:

```bash
bash run.sh --dry-run --cores 1
```

The workflow validates sample IDs, required columns, marks, assembly keys,
resolutions, paired comparisons, region viewpoints, and comparison compatibility
before scheduling jobs. Snakemake then checks file dependencies.

Also inspect:

```bash
bash run.sh --dry-run --summary
bash run.sh --list-rules
bash test.sh
```

If the dry run says a reference is missing, fix the path or run the reference
helper. Do not create an empty file to satisfy Snakemake.

## 9. Start with a pilot

The first single-run bundled sample is useful for testing wiring without launching
all 11 libraries:

```bash
bash run.sh --cores 16 results/cool/Naive_H3K27ac_rep1.mcool
```

This validates SRA download, FASTQ QC, alignment, pairtools, indexing, and cooler.
It is not a scientifically adequate single-sample analysis.

Then inspect the complete DAG again and run everything:

```bash
bash run.sh --dry-run --cores 1
bash run.sh --cores 24
```

The no-target command executes `rule all`, which is explicitly marked as the
default target. It includes matrices, peaks, primary and secondary loop calls,
QC, differential analysis, locus plots, virtual 4C, stripes, MultiQC, publication
figures, ORACLE exports, and the provenance manifest.

Snakemake is resumable. After an interruption, rerun the same command. Incomplete
outputs are rebuilt. If a stale lock remains after confirming no run is active:

```bash
bash run.sh --unlock
```

## 10. Run on SLURM

Install with `setup.sh`, then edit `profiles/slurm/config.yaml`:

- replace `compute` and `highmem` with real partition names;
- add `slurm_account` if required;
- tune memory and runtime for the site;
- place the Conda prefix on shared storage if compute nodes cannot see the default.

Run:

```bash
bash run.sh --profile profiles/slurm --jobs 100
```

The profile uses Snakemake's maintained SLURM executor, which submits and monitors
jobs through SLURM rather than a hand-written `squeue` parser.

## 11. Understand every stage

### Stage 00: SRA acquisition

`fasterq-dump` writes paired FASTQ into a per-sample temporary directory. Technical
runs are concatenated only after both mates pass integrity and record-count checks.

### Stage 01: read QC and trimming

FastQC describes base quality, adapters, composition, and duplication. fastp removes
adapters and low-quality/short reads. HiChIP often contains adapter read-through
from short molecules, so trimming is not optional.

### Stage 02: alignment and pairs

`bwa-mem2 -SP5M` uses Hi-C-aware paired-end behavior. Pairtools parses alignments,
assigns MboI fragments, coordinate-sorts, deduplicates by both ends, selects UU
contacts, removes neighbouring-fragment ligation artifacts from the contact set,
writes pre/post-filter statistics, forces BGZF compression, and creates the required
pairix `.px2` index. The
index is a declared output, so a damaged resume cannot skip it.

### Stage 03: contact matrices

Cooler loads the 5-kb unbalanced contact counts and derives all coarser resolutions.
It attempts ICE balancing at each configured resolution, then audits the HDF5
`converged` attribute. Only proven-converged weights remain in the published
`.mcool`; a nonconverged attempted weight is removed so an unaware downstream tool
cannot use it silently.

For every sample, inspect:

```text
results/qc/balance/<sample>.balance.json
results/qc/balance/<sample>.balance.tsv
```

The JSON is the machine contract and the TSV is the same resolution-by-resolution
evidence in an easy-to-read form. The states mean:

- `PASS`: the balancing weight exists and cooler reports `converged=true`;
- `WARN`: balancing ran but did not converge, so the raw matrix is retained but the
  attempted weight is not published or used;
- `NOT_ASSESSED`: the weight or convergence evidence is absent, which is unknown,
  not a pass.

A sample can therefore have, for example, `WARN` at 5 kb and `PASS` at 25 kb. This
does not invalidate its raw contacts. FitHiChIP, differential counting, HiCRep, and
stripenn use their intended raw-count inputs. Mustache is different: version 1.3.3
requires a valid `weight` column, so a non-passing balance at its configured 10-kb
resolution produces a header-only TSV and `.mustache.status.json` marked
`NOT_ASSESSED`. The balance report otherwise controls
whether each matrix-derived display or QC output uses ICE-balanced values, a labelled
raw-count fallback, or no result.

### Stage 04: anchors

The two ends of a HiChIP pair are ligation partners, not the ends of one ordinary
fragment. The workflow projects deduplicated UU alignments to individual read ends,
calls MACS3 peaks with `--nomodel --extsize 147 --keep-dup all`, and removes
blacklist overlaps. Peak calling is limited to the configured primary-assembly
chromosomes, so alternate, random, and unplaced contigs cannot leak into anchors.

Peak bases are retained when supported by at least two biological libraries in the
configured mark/tissue/protocol/enzyme assay stratum (or one when only one library
exists). Support is calculated base by base, so transitive overlaps cannot retain
unsupported flanks. The universe is shared across comparison conditions but not
across unrelated assays. The report keeps two read-end FRiP values separate: overlap
with the sample's own filtered peaks and overlap with the shared consensus anchors.
Both numerators and their denominator use the same primary autosome-plus-chrX read
population used for anchor calling. Because the assay is enriched and peaks derive
from these data, both are descriptive rather than universal pass/fail numbers.

### Stage 05: loops

FitHiChIP 11.0 models distance and 1D coverage bias in Peak-to-ALL mode. The source
archive and the extracted scripts used at runtime are checksum-verified. The
workflow reads and audits two exact FitHiChIP products:

- `interactions_FitHiC.bed`, the unthresholded, unmerged native-bin table used only
  to build the differential hypothesis universe; and
- the q-filtered `MergeNearContacts` call set used for per-sample reporting, APA,
  annotation, and visualization.

Both are normalized to explicit count/P/Q schemas using file-content compression
detection. Reporting calls remove non-primary/trans rows and blacklist overlaps,
apply the configured minimum contact count, and record every filter count. Do not
substitute the significant-call file for the all-interaction differential source.

Mustache runs independently at 10 kb as a scale-space cross-check. One-to-one caller
overlap allows the configured reciprocal-anchor tolerance (one bin by default) and
reports Jaccard plus directional support fractions. Concordance is descriptive,
never a QC gate; disagreement must be investigated, not hidden by merging callers.
The pinned caller runs through a manager-free threaded adapter so secured compute
nodes do not need Mustache's local multiprocessing socket. Output is validated,
genomically sorted, and atomically published; the status JSON records the backend.

### Stage 06: QC

The workflow produces:

- pairtools mapping, duplicate, pair-type, cis/trans, and complexity statistics;
- post-deduplication, pre-restriction-filter MboI-neighbour orientation fractions
  (not a raw-library digestion-efficiency estimate);
- per-resolution ICE convergence as JSON, TSV, and MultiQC fields;
- P(s) distance-decay curves, using a labelled raw-count fallback when needed;
- exploratory 25-kb insulation, with the same labelled fallback;
- 100-kb GC-phased E1 only when ICE balancing converged at 100 kb;
- the secondary Mustache cross-check only when ICE balancing converged at 10 kb;
- deterministic, depth-aware HiCRep SCC;
- contact-map-held-out APA against anchor/visibility/mappability-matched shifts,
  reported as a descriptive effect with a loop-bootstrap interval.

APA candidates come from sibling donors. For three-replicate groups, a candidate
normally needs both sibling call sets and must overlap sibling-only anchors. Because
those sibling FitHiChIP runs used the shared assay-stratum search space, the output
records residual cohort-level dependence instead of claiming strict independence.
The q-filtered merged reporting calls are reconciled within one native caller bin per
anchor (`apa.candidate_tolerance_bins: 1`), which tolerates donor-to-donor
representative-pixel jitter. This is intentionally different from the unthresholded
differential universe, whose hypothesis pixels must stay exact
(`differential.candidate_tolerance_bins: 0`).
Controls preserve chromosome/distance, match peak/non-peak class at both anchors
and caller-range marginal visibility, avoid blacklist bins, and require a usable
matrix window. Balanced or raw-fallback APA is `DESCRIPTIVE`, never a hard gate;
too little matched evidence is `NOT_ASSESSED`.

Cooltools cannot compute `eigs-cis` from raw counts. When 100-kb balancing is not
`PASS`, the workflow writes a header-only E1 TSV, a valid empty bigWig, and
`results/qc/compartments/<sample>.cis.eigs.status.json` with `NOT_ASSESSED`. This
preserves a stable file contract without fabricating zero E1 values.

### Stage 07: differential contacts

Each selected sample's unthresholded `interactions_FitHiC.bed` is streamed into a
normalized table. The workflow requires native 5-kb bin alignment, primary cis
chromosomes, the exact configured 20-kb to 3-Mb range, and no blacklist overlap. It
then constructs a condition-blind exact-pixel universe: a pixel is eligible only
when its raw contact count is at least `min_count` in at least `min_samples`
libraries. The universe does not inspect q-values or condition labels and does not
merge, widen, or tolerate neighbouring pixels.

The same exact pixels are recounted from every sample's blacklist-filtered
unbalanced cooler. The workflow validates unique sample/pixel keys, complete sample
coverage, and equality between a source FitHiChIP count and its cooler recount.
Missing or duplicated rows fail; they are never filled with an invented zero.
`hypothesis_universe.json` records source files, thresholds, selected samples,
per-sample source counts, and the zero-tolerance contract.

The bundled paired model is:

```text
counts ~ subject_id + condition
```

pyDESeq2 estimates size factors and dispersion, tests the condition coefficient,
and controls FDR. The output includes the hypothesis count, low-count filter,
metadata, fitted factors, exact design formula, full Wald results, and
`paired_effects.tsv` so donor directions can be inspected. `design.json` reports
`n_complete_pairs`, `publication_eligible`, and either `STANDARD_INFERENCE` or
`PILOT_UNDERPOWERED`. The bundled two-pair contrasts are pilots even if an adjusted
p-value happens to cross a conventional threshold.

### Stage 08: locus visualization

Each composite includes a balance-aware contact heatmap, exploratory insulation,
peaks, filtered loop arcs, GENCODE gene models, coordinates, and resolution.
Heatmaps use ICE values only when that resolution is `PASS`; otherwise they use raw
counts and say so in the plot. Temporary sidecars are region-specific, avoiding
races when many regions render together.

Virtual 4C reads one matrix row from the explicit viewpoint, avoiding a dense
whole-chromosome allocation. Its title and y-axis state `ICE-balanced` or
`raw-count fallback` together with the balance status. It is a visual profile, not
a differential test.

### Stage 09: ORACLE export

Each sample is exported as PyTorch Geometric and HDF5 representations with multiple
bin scales connected by explicit `contained_by`/`contains` hierarchy edges. Loop
and genomic-adjacency relations remain separate. FitHiChIP count, raw P value, and
adjusted Q/FDR remain distinct, and coincident coarse loop edges are consolidated
with their fine-call multiplicity recorded. Coarse
insulation/E1 intervals are overlap-projected across all covered fine bins; missing
values are represented explicitly before the documented export fill policy.

The node-channel order stays fixed even when a structural track is unavailable.
Read `node_feature_availability` in the per-sample manifest before using `x`:
insulation records whether it came from ICE or the raw-count fallback, while E1 is
`available: false` when 100-kb balancing was not `PASS`. Every false
`x_observed_mask` cell, including all channels at a blacklisted node, is zero-filled;
those zeros are not measured biological values. The manifest also
records and hashes the balance report and the E1 status evidence.

Treat these node features as a structural prototype. A future multimodal integration
should replace or augment peak-overlap indicators with continuous measured signal.

### Stages 10–12: report, figures, provenance

Stable flat JSON remains the machine contract. Lightweight `_mqc.json` companions
render HiChIP-specific sample, APA, and differential-design sections in MultiQC
without rerunning analyses. The sample table includes the overall balance state and
converged, nonconverged, or missing resolutions, so a `WARN` is visible without
opening the HDF5 file. The HTML is self-contained.

Five cohort figures are written as vector PDF and 400-dpi PNG:

1. library/contact-map QC;
2. same-mark reproducibility context and HiCRep;
3. primary-role same-mark loop yield, depth, span, and held-out APA;
4. paired differential evidence, including donor-level directions and Wald
   uncertainty with an explicit pilot/publication status;
5. stripes for primary-role libraries, with configured demonstration exclusions
   stated on the figure and retained in QC/raw tables.

Figure inputs are required schema contracts: a missing column or malformed required
table stops rendering instead of producing a plausible no-data panel. Empty but
valid results remain visibly labelled.

The final run manifest embeds and hashes the effective configuration (including
command-line overrides), actual reference files, resolved packages for portable
Conda rule environments, QC tables, designs, report/figure files, and ORACLE
manifests, and records the Git revision.

## 12. Read QC results correctly

Open:

```bash
# Linux desktop
xdg-open results/multiqc/multiqc_report.html

# or copy the self-contained HTML to your workstation
```

### Balance convergence

Start with the balance fields in MultiQC's sample table, then open the per-sample
JSON or TSV when any resolution is not `PASS`. Interpret the states literally:

- `PASS` is positive evidence that the requested ICE iteration converged.
- `WARN` says only that ICE did not converge at one or more resolutions. It does not
  say that the raw contact map or raw-count loop analysis failed.
- `NOT_ASSESSED` means the convergence result is unavailable and must remain unknown.

Next, check the `normalization` and `balance_status` columns or plot labels for each
matrix-derived output. Expected-cis, insulation, locus heatmaps, and virtual 4C can
fall back to raw counts. These remain useful diagnostics but are more sensitive to
coverage and mappability and are not equivalent to ICE-balanced tracks. APA is
descriptive under either normalization and preserves that label; E1 is
`NOT_ASSESSED` because no scientifically valid raw-count fallback exists.

### Valid-pair yield

```text
deduplicated valid-ligation UU contacts / all raw sequenced read pairs
```

The denominator comes from fastp's before-filtering read count divided by two for
paired-end data. The report also gives a descriptive post-trim denominator. This is
not UU / UU, which would always be 100%.

### Duplicate rate

```text
duplicate mapped pairs / all mapped pairs
```

It uses pairtools' mapped-pair denominator. A high rate can indicate limited library
complexity or over-sequencing. Interpret it alongside absolute usable contacts.

### Cis fraction

```text
cis deduplicated valid-ligation UU contacts / all deduplicated valid-ligation UU contacts
```

The shipped 0.70 threshold is a project screen, not a law across all cell types,
protocols, or depths.

### Loop count

Loop yield depends on depth, mark, anchor number, and biology. Never interpret a
bar chart without the depth plot and replicate agreement. A low-depth sample with
few loops is not evidence of absent biology.

### HiCRep

Only compare libraries in the same configured biological condition, tissue, mark,
and protocol group. A pair below the configured contact floor is depth-confounded
and cannot decide reproducibility. For the remaining pairs, a
sample and its full replicate group are `PASS` only when every relevant pair clears
the SCC threshold, `FAIL` when every pair is below it, `DISCORDANT` when results are
mixed, and `NOT_ASSESSED` when none qualify. Minimum, mean, best, and group-median
SCC values remain useful descriptions, but none is the gate. Stable downsampling
seeds make reruns exact.

### APA

- `DESCRIPTIVE`: a matched-control effect and loop-bootstrap interval were computed;
- `NOT_ASSESSED`: too few sibling-supported loops, usable windows, or matched controls.

The interval describes within-dataset loop-resampling uncertainty. It is not a
p-value or a universally calibrated library-quality threshold.

## 13. Interpret differential output

For each comparison, inspect:

```text
results/loops/<sample>/<sample>.interactions_FitHiC.all.audit.json
results/diff/<comparison>/hypothesis_universe.json
results/diff/<comparison>/candidate_support.tsv
results/diff/<comparison>/design.json
results/diff/<comparison>/paired_effects.tsv
results/diff/<comparison>/differential_loops.tsv
results/diff/<comparison>/ma_plot.png
results/diff/<comparison>/volcano.png
```

Before reading significant rows, verify:

1. the design JSON contains the intended `include_subjects` subset,
   `~ subject_id + condition`, and an honestly interpreted `analysis_status`;
2. all pairs are biological, not technical, replicates;
3. no arm is systematically shallower or lower-complexity;
4. the universe manifest says `fithichip_all_interactions`, tolerance zero, and the
   expected `min_count`/`min_samples` abundance filter;
5. all-interaction audits and strict recount checks cover the expected samples and
   exact native pixels;
6. effect directions are supported by per-sample counts and `paired_effects.tsv`;
7. important loci are visible in the contact map and sibling-donor evidence;
8. language acknowledges that occupancy and 3D contact are not separated.

With only two primary pairs, `analysis_status` is `PILOT_UNDERPOWERED`. An exhaustive
paired label-swap test has four assignments and cannot provide a conventional 0.05
empirical p-value. Do not promote Wald threshold hits to ordinary significant
discoveries or claim that a tiny permutation space “validates” the result. Prefer
additional donors and orthogonal validation.

## 14. Troubleshooting

### `setup.sh --check` reports fewer than 19 environments

Run:

```bash
bash setup.sh
```

Do not hand-edit `.snakemake/conda`.

### A reference is missing

```bash
bash prepare_references.sh hg38
```

If using custom files, update every related path in `config/genome.yaml` and keep
FASTA, annotation, blacklist, chromosome sizes, indexes, and digest on one assembly.

### SRA download stops

Rerun the same command. The rule removes an incomplete one-mate cache and only
publishes paired, structurally valid gzip files. Check network, quota, and temporary
storage in `data/sra_tmp/`.

### Pairix index is missing

The `.px2` file is a declared output. Rerun the target; Snakemake will rebuild the
deduplicated pairs rule rather than allowing cooler to fail later.

### FitHiChIP cannot download on a compute node

Run `bash setup.sh` on the login node first. It stores the checksum-verified source
under `.cache/downloads/`, which is ignored by Git but available to the rule.

### MultiQC lacks loop or APA sections

Confirm files ending `_loop_qc_mqc.json` and `_apa_mqc.json` exist under `results/qc`
and the differential MultiQC companions exist under `results/diff`:

```bash
bash run.sh --cores 4 multiqc --force
```

The stable flat JSON alone is intentionally not parsed as custom content.

### Balance status is `WARN` or `NOT_ASSESSED`

Open `results/qc/balance/<sample>.balance.tsv` and identify the affected resolution.
Check usable cis depth, sparsity, blacklist/mappability patterns, and the cooler log
before changing balancing parameters. A very fine HiChIP matrix may legitimately be
too sparse while coarser resolutions pass. Do not copy weights from another sample,
restore a removed nonconverged weight, or label a raw fallback as balanced. If the
resolution is essential, increase usable depth or justify a revised, cohort-wide
balancing policy and rerun every comparable sample.

### A plotting rule fails for one empty loop set

A schema-valid empty biological result renders an explicit no-data panel. A malformed
required BEDPE/table header, missing sample, or assembly mismatch must fail. Fix the
input contract; do not add `|| true` or a broad exception to publish a misleading
figure.

### A run was killed

Check the Snakemake and per-rule logs, increase the relevant profile resource, and
rerun. Snakemake will reuse complete outputs. Do not delete all results as a first
response.

## 15. Final review checklist

Before sharing results, confirm:

- [ ] `bash setup.sh --check` passes.
- [ ] `bash run.sh --dry-run --cores 1` reports nothing to do.
- [ ] `bash test.sh` passes.
- [ ] Every biological group has the intended independent donors.
- [ ] Valid yield, duplicate rate, cis fraction, both FRiP definitions, depth, and
      post-deduplication restriction QC were inspected using their stated
      denominators.
- [ ] Balance JSON/TSV and MultiQC were reviewed for every required resolution.
- [ ] Every raw-count fallback is labelled, APA remains descriptive, and a
      non-passing 100-kb balance leaves E1 `NOT_ASSESSED`.
- [ ] A non-passing Mustache resolution is `NOT_ASSESSED`; no invalid weight was recreated.
- [ ] HiCRep depth-confounded pairs were not treated as failures or successes.
- [ ] APA used contact-map-held-out sibling candidates, sibling-only anchors,
      matched controls, and reports residual search-space dependence.
- [ ] Differential `hypothesis_universe.json` proves unthresholded all-interaction
      input, native pixels, zero tolerance, and the configured abundance filter.
- [ ] Differential `design.json` contains the intended pairing, factors, complete-pair
      count, and honestly interpreted `analysis_status`; paired effects were reviewed.
- [ ] Loop/stripe differences were checked against depth and mark.
- [ ] Locus figures include gene models and explicit viewpoints.
- [ ] E1/insulation language follows `docs/INTERPRETATION.md`.
- [ ] ORACLE consumers check `node_feature_availability`, `x_observed_mask`, and
      `blacklist_mask` rather than treating filled zeros as measurements.
- [ ] `results/provenance/run_manifest.json` exists and records actual reference
      hashes plus resolved portable rule-environment packages.
- [ ] No FASTQ, BAM, matrices, results, caches, notes, worktrees, or local analysis
      artifacts are staged for Git.

At that point the workflow is complete. The scientific work continues with
orthogonal validation, careful effect interpretation, and transparent reporting of
the experimental design and limitations.
