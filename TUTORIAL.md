# ORACLE HiChIP Pipeline — Step-by-Step Tutorial

This tutorial walks you through the complete pipeline from raw FASTQ files to
ORACLE-ready graph objects, explaining **what each step does**, **why it matters**,
and **what to look for** in the output.

---

## Table of Contents

1. [What is HiChIP?](#1-what-is-hichip)
2. [Install the environment](#2-install-the-environment)
3. [Download reference files](#3-download-reference-files)
4. [Prepare your sample sheet](#4-prepare-your-sample-sheet)
5. [Configure the pipeline](#5-configure-the-pipeline)
6. [Run the pipeline](#6-run-the-pipeline)
7. [Step-by-step walkthrough](#7-step-by-step-walkthrough)
   - [Step 01 — Raw read QC](#step-01--raw-read-qc)
   - [Step 02 — Alignment & pair extraction](#step-02--alignment--pair-extraction)
   - [Step 03 — Contact matrix](#step-03--contact-matrix)
   - [Step 04 — 1D peak calling](#step-04--1d-peak-calling)
   - [Step 05 — Loop calling (FitHiChIP)](#step-05--loop-calling-fithichip)
   - [Step 06 — QC the loops](#step-06--qc-the-loops)
   - [Step 07 — Differential loops](#step-07--differential-loops-optional)
   - [Step 08 — Visualisation](#step-08--visualisation)
   - [Step 09 — ORACLE export](#step-09--oracle-export)
   - [Step 10 — MultiQC report](#step-10--multiqc-report)
8. [Understanding the outputs](#8-understanding-the-outputs)
9. [Common problems and fixes](#9-common-problems-and-fixes)
10. [Frequently asked questions](#10-frequently-asked-questions)

---

## 1. What is HiChIP?

HiChIP combines **Hi-C** (genome-wide chromatin contacts) with **ChIP-seq**
(protein-of-interest immunoprecipitation). The result is a sparse contact map
**enriched for interactions anchored at the immunoprecipitated protein** (e.g.
H3K27ac for active enhancers, CTCF for structural loops).

```
Library prep overview
─────────────────────
Cells → crosslink → restrict (MboI/DpnII) → ligate → immunoprecipitate
     → sonicate → paired-end sequence → YOUR FASTQ FILES
```

**What the pipeline produces:**

| Output | Format | Used for |
|--------|--------|---------|
| Trimmed FASTQs | `.fastq.gz` | Alignment input |
| Contact pairs | `.pairs.gz` | Matrix + loop calling |
| Contact matrix | `.mcool` | Visualisation, ORACLE |
| 1D peaks | `.bed` | Loop anchors |
| Loops | `.bed` (BEDPE-style) | Differential analysis, ORACLE |
| QC report | `.html` (MultiQC) | Library quality assessment |
| ORACLE graph | `.pt` / `.h5` | Foundation model input |

---

## 2. Install the environment

### Prerequisites

- Linux (Ubuntu 20.04+ or equivalent; macOS not supported for all tools)
- ≥ 64 GB RAM (bwa-mem2 index build needs ~70 GB for hg38)
- ≥ 200 GB disk for references + results

### One-command setup

```bash
# Clone the repository
git clone https://github.com/challagandla/oracle-hichip.git
cd oracle-hichip

# Build the conda environment (≈ 10–20 min)
bash setup_env.sh

# Activate
mamba activate oracle-hichip

# Confirm key tools are available
snakemake --version          # should print 8.x
cooler --version             # should print 0.10.x
pairtools --version          # should print 1.1.x
```

> **GPU note**: The environment includes PyTorch + CUDA 12.1 for the ORACLE export step.
> If you are on a CPU-only machine, remove the `pytorch-cuda` line from `environment.yml`
> before running `setup_env.sh`.

---

## 3. Download reference files

Before running the pipeline you need a reference genome, aligner index,
chromosome sizes, blacklist, GTF, and restriction digest BED.
A download script handles all of this:

```bash
# Download hg38 references (takes ~60 min + bwa-mem2 indexing time)
bash resources/download_resources.sh hg38

# Download multiple assemblies
bash resources/download_resources.sh hg38 mm10
```

The script will tell you the exact paths to put in `config/genome.yaml`.

**What gets downloaded:**

| File | Source |
|------|--------|
| `GRCh38.primary_assembly.genome.fa` | GENCODE / EBI |
| `gencode.v46.primary_assembly.annotation.gtf.gz` | GENCODE |
| `hg38-blacklist.v2.bed.gz` | Boyle Lab |
| `hg38.chrom.sizes` | UCSC |
| `bwamem2_index/` | Built locally from FASTA |
| `bwa_index/` | Built locally from FASTA |
| `MboI.digest.hg38.bed.gz` | Generated locally from FASTA |

---

## 4. Prepare your sample sheet

Edit `config/samples.tsv`. Each row is one HiChIP sample (one library).

```tsv
sample_id         subject_id  tissue   disease           replicate  mark     fastq_r1                           fastq_r2                           batch      library_protocol  notes
OV_tumor_K27ac_R1 OV001       tumor    ovarian_HGSOC     1          H3K27ac  data/fastq/OV001_K27ac_R1_R1.fastq.gz  data/fastq/OV001_K27ac_R1_R2.fastq.gz  2026Q2_A   HiChIP_v2         Pilot cohort
OV_tumor_K27ac_R2 OV001       tumor    ovarian_HGSOC     2          H3K27ac  data/fastq/OV001_K27ac_R2_R1.fastq.gz  data/fastq/OV001_K27ac_R2_R2.fastq.gz  2026Q2_A   HiChIP_v2         Pilot cohort
```

**Column reference:**

| Column | Required | Description |
|--------|----------|-------------|
| `sample_id` | ✅ | Unique identifier; used in all output filenames |
| `subject_id` | ✅ | Donor/patient ID; used for replicate grouping in HiCRep |
| `tissue` | ✅ | Tissue type (e.g. `tumor`, `pbmc`, `cortex`) |
| `disease` | ✅ | Disease label (or `healthy`) |
| `replicate` | ✅ | Integer replicate number within subject + mark |
| `mark` | ✅ | Histone mark or factor (e.g. `H3K27ac`, `CTCF`). Must match a key in `config.yaml:macs3.marks` |
| `fastq_r1` | ✅ | Path to R1 FASTQ (absolute or relative to repo root) |
| `fastq_r2` | ✅ | Path to R2 FASTQ |
| `batch` | ✅ | Sequencing batch label; used to detect batch effects in differential analysis |
| `library_protocol` | optional | `HiChIP_v2`, `CUT_Tag_HiChIP`, etc. |
| `notes` | optional | Free text |

> ⚠️ **Adapter sequences** — if you mix `HiChIP_v2` and `CUT_Tag_HiChIP` in the same
> run, remember that CUT&Tag uses Tn5 mosaic adapters, not TruSeq.
> Edit `config.yaml:fastp.adapter_sequence` accordingly (or split into separate pipeline runs).

---

## 5. Configure the pipeline

Open `config/config.yaml` and check these settings before your first run:

```yaml
# 1. Point to your sample sheet
samples_tsv: config/samples.tsv

# 2. Choose the genome assembly (must match a key in config/genome.yaml)
assembly: "hg38"    # or T2T-CHM13 / mm10

# 3. Set thread counts to match your machine
threads:
  bwa: 24           # set to the number of cores on your alignment node
  fithichip: 16

# 4. Check MACS3 genome size
macs3:
  genome_size: "hs"   # hs for human, mm for mouse

# 5. Leave differential.comparisons empty until you have matched controls
differential:
  comparisons: []
```

For most users the defaults are good. The most common changes are thread counts
and the `differential.comparisons` block once you have case/control pairs.

---

## 6. Run the pipeline

```bash
mamba activate oracle-hichip

# Always do a dry-run first — see exactly what will run
snakemake -s workflow/Snakefile -n --configfile config/config.yaml

# Visualise the DAG (requires graphviz)
snakemake -s workflow/Snakefile --dag --configfile config/config.yaml | dot -Tsvg > dag.svg

# Full run on a local workstation (32 cores)
snakemake -s workflow/Snakefile --cores 32 --configfile config/config.yaml --use-conda

# Run only QC steps (useful for initial library assessment)
snakemake -s workflow/Snakefile --cores 16 qc_raw --configfile config/config.yaml

# Run everything up to and including loop calling
snakemake -s workflow/Snakefile --cores 32 loops_fithichip --configfile config/config.yaml

# On a SLURM cluster
snakemake -s workflow/Snakefile --profile profiles/slurm --configfile config/config.yaml

# Resume after a failed run (Snakemake automatically reruns only failed/missing outputs)
snakemake -s workflow/Snakefile --cores 32 --configfile config/config.yaml --rerun-incomplete
```

**Convenience targets** (run only part of the pipeline):

| Target | What it builds |
|--------|---------------|
| `qc_raw` | FastQC + fastp reports |
| `align_pairs` | Trimmed FASTQs → `.pairs.gz` |
| `cool_matrix` | `.mcool` contact matrices |
| `peaks` | MACS3 1D peak BEDs |
| `loops_fithichip` | FitHiChIP loop BEDs |
| `loop_qc` | Per-sample QC JSON summaries |
| `export_oracle` | ORACLE `.pt` + `.h5` graph files |
| `multiqc` | Aggregated HTML QC report |
| `differential_all` | DE loops for all configured comparisons |

---

## 7. Step-by-step walkthrough

### Step 01 — Raw read QC

**What it does**: FastQC checks read quality; fastp trims adapter sequences
and low-quality bases.

**Why it matters**: HiChIP libraries frequently contain adapter readthrough
(the insert is shorter than the read length, so you sequence into the adapter).
If you skip trimming, pairtools will misclassify adapter-contaminated reads as
chimeric pairs and inflate your junk-pair fraction.

**Outputs**:
```
results/qc/fastqc_raw/<sample>_R1_fastqc.html   ← open in browser
results/qc/fastp/<sample>.fastp.html             ← open in browser
results/trimmed/<sample>_R1.trim.fastq.gz        ← used in Step 02
```

**What to check**:
- Per-base quality should be ≥ Q20 for ≥ 70% of reads after trimming
- Adapter content column in FastQC should show low or no adapter after trimming
- fastp HTML shows the adapter detection and trimming summary

---

### Step 02 — Alignment & pair extraction

**What it does**:
1. Aligns trimmed reads with **bwa-mem2 -SP5M** (canonical Hi-C mode)
2. Parses the BAM with **pairtools** to extract contact pairs
3. Sorts, deduplicates pairs, keeps only `UU` (unique-unique) pairs
4. Indexes the `.pairs.gz` with pairix

**Why `-SP5M`?**
- `-S` skip mate rescue (prevents false-positive chimeric pairs)
- `-P` skip pairing (pairtools handles this)
- `-5` report secondary alignments
- `-M` mark short split alignments as secondary

**Why pairtools, not HiC-Pro?**
pairtools is actively maintained, emits the standard `.pairs.gz` format, and
defines duplicates correctly — as pairs where both ends match, not individual reads.

**Outputs**:
```
results/pairs/<sample>.dedup.pairs.gz      ← clean contact pairs
results/pairs/<sample>.dedup.pairsam.gz   ← includes SAM for peak calling
results/qc/pairtools/<sample>.dedup.stats.txt
```

**What to check**:
- `pairtools stats` output: `total_nodups` / `total` = valid pair yield; target ≥ 25%
- Duplicate fraction: target < 50%; > 50% = undersampled library, consider resequencing
- `UU` fraction of total: target > 60%

---

### Step 03 — Contact matrix

**What it does**:
1. Bins pairs at the finest resolution (5 kb) with **cooler cload pairix**
2. Zoomifies to all configured resolutions (5 kb → 2.5 Mb) with **cooler zoomify**
3. ICE-balances every resolution during zoomify

**Why `.mcool` and not `.hic`?**
`.mcool` (HDF5-backed, multi-resolution) is the modern standard used by cooltools,
HiGlass, and ORACLE. It is more efficient to query than `.hic`. If you need `.hic`
for Juicebox, use `cooler dump` + Juicer tools to convert.

**Outputs**:
```
results/cool/<sample>.base.cool   ← 5 kb single-resolution
results/cool/<sample>.mcool       ← all resolutions, balanced
```

**What to check**:
- Matrix should be non-empty: `cooler info results/cool/<sample>.mcool`
- Cis/trans ratio from pairtools stats: cis ≥ 70% of valid pairs

---

### Step 04 — 1D peak calling

**What it does**:
1. Splits the deduped pairsam back to a 1D BAM (both read ends)
2. Runs **MACS3** to call peaks in narrow mode (H3K27ac, H3K4me3, CTCF)
   or broad mode (H3K27me3, H3K36me2, H3K36me3)

**Why call peaks from the HiChIP reads themselves?**
The immunoprecipitation step in HiChIP enriches for the protein of interest.
Using a separate ChIP-seq aliquot is sub-optimal because it comes from a different
ligation pool. The 1D read distribution from HiChIP directly reflects where your
antibody pulled down.

**Outputs**:
```
results/peaks/<sample>_peaks.bed   ← clean BED3 used by FitHiChIP
results/peaks/raw/<sample>_peaks_macs.done
```

**What to check**:
- Number of peaks: typically 30,000–150,000 for H3K27ac; fewer for CTCF
- Open the narrowPeak/broadPeak in a genome browser and verify it makes biological sense

---

### Step 05 — Loop calling (FitHiChIP)

**What it does**:
1. Converts the `.pairs.gz` to FitHiChIP's `allValidPairs` format
2. Generates a per-sample FitHiChIP config
3. Runs **FitHiChIP** (Peak-to-ALL mode, 5 kb bins, 20 kb–3 Mb, FDR < 0.01)
4. Cross-checks with **mustache** as a sanity check

**Why FitHiChIP and not a generic loop caller?**
FitHiChIP models the **1D ChIP enrichment bias** using spline regression. HiChIP
loops are called in a biased contact map (peaks are more connected simply because
they are enriched in reads), so without bias correction you get thousands of
false-positive "loops" that are just ChIP peaks near each other.

**Outputs**:
```
results/loops/<sample>/<sample>.interactions_FitHiC_Q0.01.bed   ← primary loops
results/loops/<sample>/<sample>.mustache.tsv                     ← cross-check
```

**What to check**:
- Loop count: ≥ 1,000 loops at q < 0.01 per sample is the QC floor
- Fewer loops → poor library, wrong antibody, or overly strict thresholds
- Cross-check: loops called by both FitHiChIP and mustache are the most reliable

---

### Step 06 — QC the loops

**What it does**: Aggregates all QC metrics into per-sample JSON summaries:
- **P(s) distance decay** — `cooltools expected-cis`
- **Insulation scores** — TAD boundary scores at 25 kb
- **A/B compartments** — first eigenvector at 100 kb
- **APA** — Aggregate Peak Analysis; pile-up of contact maps at loop positions
- **HiCRep** — stratum-adjusted correlation between biological replicates

**APA explained**:

```
  APA score = contact_at_loop_centre / mean(corner_contacts)
  ≥ 1.5 vs random-shift controls = PASS
```

If the score is < 1.5 it means your loops don't show contact enrichment above
background, which indicates either poor loop calling or a poor-quality library.

**HiCRep SCC** — read this one carefully, it is not what it looks like:
- Replicates are grouped on **cell type + mark**, not on donor. Several conditions
  taken from the same donors would otherwise land in one "replicate" group.
- SCC is **depth-dominated**: HiCRep downsamples the deeper library to the shallower,
  so the pair is only as good as its worst member. On GSE101498 at 25 kb, two
  libraries from *different* cell types (90.8M and 47.6M contacts) score **0.743**,
  while two genuine replicates score **0.221** because one holds 3.0M.
- So it is a sanity floor (`hicrep.threshold_pass`), not proof that replicates are
  tighter than conditions. Any pair whose shallower member is below
  `hicrep.min_contacts_for_scc` is reported as `depth_confounded` and does not decide
  PASS/FAIL.
- Single-replicate samples are reported as `NOT_ASSESSED` (not PASS)

**Outputs**:
```
results/qc/expected/<sample>.expected.cis.tsv
results/qc/insulation/<sample>.insulation.tsv
results/qc/compartments/<sample>.cis.eigs.tsv
results/qc/compartments/<sample>.E1.bw         ← bigWig for browser
results/qc/apa/<sample>.apa.png
results/qc/apa/<sample>.apa.json
results/qc/hicrep/<sample>.hicrep.json
results/qc/loop_qc/<sample>.json               ← overall pass/fail
results/qc/loop_qc/<sample>.md
```

**What to check — QC thresholds**:

| Metric | PASS threshold |
|--------|---------------|
| Valid pair yield | ≥ 25% |
| Duplicate fraction | ≤ 50% |
| Cis fraction | ≥ 70% |
| Loop count | ≥ 1,000 at q < 0.01 |
| APA score | ≥ 1.5 vs random-shift controls (distance-matched; not the window corners) |
| HiCRep SCC | ≥ `hicrep.threshold_pass`, and only for pairs above `hicrep.min_contacts_for_scc` — below that it grades depth, not concordance |

---

### Step 07 — Differential loops (optional)

**What it does**: For explicitly configured case/control groups,
counts valid pairs supporting each loop in each sample, then runs
**pyDESeq2** (pseudobulk-style) on the loop-by-sample count matrix.

**How to enable it**: Add a comparison block in `config/config.yaml`:
```yaml
differential:
  fdr: 0.05
  log2fc_min: 1.0
  method: "pyDESeq2"
  comparisons:
    - name: ovarian_tumor_vs_adjacent_K27ac
      mark: H3K27ac
      case_filter:
        tissue: tumor
        disease: ovarian_HGSOC
        library_protocol: HiChIP_v2
      control_filter:
        tissue: adjacent_normal
        disease: non_cancer
        library_protocol: HiChIP_v2
```

> ⚠️ **Never compare samples with different marks, tissues, or protocols.**
> The pipeline will raise an error if case and control groups are empty, but it
> cannot detect semantic mismatches (e.g., H3K27ac vs H3K36me2).

**Requirements**:
- ≥ 3 biological replicates per group (more is better)
- Replicates must share the same `mark`, `tissue`, and `library_protocol`
- Each replicate must be a separate row in `samples.tsv`

**Outputs**:
```
results/diff/<comparison>/differential_loops.tsv   ← log2FC, padj, coordinates
results/diff/<comparison>/volcano.png
results/diff/<comparison>/ma_plot.png
results/diff/<comparison>/design.json
```

---

### Step 08 — Visualisation

**What it does**: Generates publication-grade static figures for configured
regions of interest using **pyGenomeTracks**:
- HiChIP contact heatmap (log1p scale)
- Insulation score track
- 1D peak track (mark-specific)
- FitHiChIP loop arcs

Also generates:
- **Virtual 4C**: contact profile from a viewpoint (anchor) of interest
- APA aggregate plots (from Step 06)

**Configure regions** in `config/config.yaml`:
```yaml
viz:
  regions:
    - { name: "MYC_locus",   chrom: "chr8",  start: 126500000, end: 129500000 }
    - { name: "MYCN_locus",  chrom: "chr2",  start: 16000000,  end: 17000000 }
    - { name: "BCL2_locus",  chrom: "chr18", start: 60500000,  end: 63500000 }
```

**Outputs**:
```
results/viz/<sample>_<region>.png               ← composite figure
results/viz/virtual_4c/<sample>_<region>.v4c.bw
results/viz/virtual_4c/<sample>_<region>.v4c.png
```

---

### Step 09 — ORACLE export

**What it does**: Converts each sample into the ORACLE Chromatin Operating
System (COS) graph format for foundation-model training.

**Graph structure**:
```
Nodes  = genomic bins at 5 kb / 25 kb / 100 kb / 1 Mb
Edges  = FitHiChIP loop contacts + genomic adjacency edges
Node features (current prototype):
  [0] peak_overlap_count_per_kb  — MACS3 peak overlap, normalised by bin size
  [1] insulation                 — cooltools insulation score
  [2] E1_eigenvector             — A/B compartment eigenvector
Edge features:
  [0] loop_score                 — FitHiChIP confidence score
  [1] loop_fdr                   — adjusted p-value
  [2] genomic_distance_bp        — linear distance between anchors
```

> ⚠️ **`peak_overlap_count_per_kb` is a prototype feature**.
> It counts how many MACS3 peak intervals overlap each bin — it is NOT continuous
> per-mark ChIP/CUT&Tag signal. Full ORACLE COS will merge continuous bigWig tracks
> from sister ATAC/CUT&Tag/RNA pipelines.

**Outputs**:
```
results/oracle_cos/<sample>.pt             ← PyTorch Geometric HeteroData
results/oracle_cos/<sample>.h5             ← HDF5 mirror (no PyTorch needed)
results/oracle_cos/<sample>.manifest.json  ← provenance: SHA-256 of inputs,
                                              feature names, caveats
```

**Load the graph in Python**:
```python
import torch
data = torch.load("results/oracle_cos/MY_SAMPLE.pt")

# Multi-resolution node features at 25 kb
x_25kb = data["bin_res_25000"].x          # shape: (n_bins, 3)

# Loop edges at 25 kb
ei = data[("bin_res_25000", "contact", "bin_res_25000")].edge_index
ea = data[("bin_res_25000", "contact", "bin_res_25000")].edge_attr
ek = data[("bin_res_25000", "contact", "bin_res_25000")].edge_kind
# edge_kind: 1 = loop, 0 = genomic adjacency
```

---

### Step 10 — MultiQC report

**What it does**: Aggregates all FastQC, fastp, pairtools, loop QC, APA,
and HiCRep metrics into a single interactive HTML report.

```bash
# Open in browser
xdg-open results/multiqc/multiqc_report.html
```

The report includes:
- Read QC table (FastQC pass/fail flags per sample)
- Adapter trimming statistics (fastp)
- Alignment + pairs statistics (pairtools)
- **Loop QC summary table** (valid pair yield, duplicate %, cis fraction,
  loop count, APA score, HiCRep SCC, overall PASS/FAIL status)

---

## 8. Understanding the outputs

### Output directory structure
```
results/
├── trimmed/              ← fastp-trimmed FASTQs
├── qc/
│   ├── fastqc_raw/       ← FastQC HTMLs
│   ├── fastp/            ← fastp JSON + HTML
│   ├── pairtools/        ← pairs statistics
│   ├── expected/         ← P(s) distance decay
│   ├── insulation/       ← TAD boundary scores
│   ├── compartments/     ← A/B eigenvectors + bigWigs
│   ├── apa/              ← APA plots + score JSON
│   ├── hicrep/           ← replicate concordance JSON
│   └── loop_qc/          ← aggregated QC JSON + Markdown per sample
├── pairs/                ← pairtools .pairs.gz + pairsam
├── bam_1d/               ← 1D BAM for MACS3
├── cool/                 ← .base.cool and .mcool matrices
├── peaks/                ← MACS3 peak BEDs
├── loops/                ← FitHiChIP + mustache loop BEDs
├── stripes/              ← stripenn architectural stripes + summary
├── diff/                 ← differential loop results
├── viz/                  ← pyGenomeTracks figures + virtual 4C
├── figures/              ← cohort publication figures (PDF + PNG)
├── oracle_cos/           ← ORACLE .pt / .h5 / manifest.json
└── multiqc/              ← multiqc_report.html
```

### Key files to archive

After a successful run, the minimal set to archive (excluding large intermediate files):

```bash
results/cool/*.mcool              # Contact matrices
results/peaks/*_peaks.bed         # 1D peaks
results/loops/*/*interactions*.bed # Called loops
results/qc/loop_qc/*.json         # QC summaries
results/oracle_cos/*.pt           # ORACLE graphs
results/oracle_cos/*.manifest.json # Provenance
results/multiqc/multiqc_report.html
config/                           # All configs
environment.lock.yml              # Exact environment
```

---

## 9. Common problems and fixes

### "fithichip: command not found" / "FitHiChIP_HiCPro.sh: command not found"

The pipeline automatically tries both `fithichip` (bioconda) and
`FitHiChIP_HiCPro.sh` (legacy). If neither is found:

```bash
mamba install -c bioconda fithichip
# or, for the legacy version:
conda install -c bioconda fithichip=1.0.*
```

### "pairix: pairs file not indexed"

The dedup rule should index the `.pairs.gz` automatically with `pairix -f -p pairs`.
If it fails (e.g. pairs file is corrupt):

```bash
pairix -f -p pairs results/pairs/<sample>.dedup.pairs.gz
```

### Zero loops from FitHiChIP

Common causes:
1. **Peak file is empty** — check `results/peaks/<sample>_peaks.bed`; if empty, MACS3 failed silently.
2. **Too few valid pairs** — check pairtools stats; if valid pair yield < 10%, the library needs resequencing.
3. **Wrong `bin_size`** — at very low sequencing depth (< 100M total reads), 5 kb resolution is too fine. Try 10 kb: `fithichip.bin_size: 10000`.

### MultiQC report is missing loop QC tables

Ensure you pass `--config config/multiqc_config.yaml` (this is now handled automatically
by the pipeline, but if you run MultiQC manually, include it):

```bash
multiqc --config config/multiqc_config.yaml -f -o results/multiqc results/
```

### bwa-mem2 index fails ("out of memory")

bwa-mem2 indexing needs ~70 GB RAM for hg38. If your machine doesn't have enough,
either use the legacy bwa index:

```yaml
# config/config.yaml
bwa:
  use_bwamem2: false
```

Or run index building on a high-memory node:

```bash
bwa-mem2 index -p resources/hg38/bwamem2_index/GRCh38.fa resources/hg38/GRCh38.fa
```

### Snakemake says "missing output" after apparent success

Usually caused by a `|| true` suppressing an error in older versions.
The current codebase removes all `|| true` from output-producing steps.
Run with `--rerun-incomplete` to force rerun:

```bash
snakemake -s workflow/Snakefile --cores 32 --configfile config/config.yaml --rerun-incomplete
```

---

## 10. Frequently asked questions

**Q: Can I use this pipeline with data from Arima HiChIP, Dovetail, or Phase Genomics kits?**

Yes. All kits use MboI or a compatible enzyme. The key difference is the adapter
sequence in `config/config.yaml:fastp.adapter_sequence`. Arima uses standard Illumina
TruSeq adapters; Dovetail and Phase Genomics also use TruSeq by default (verify in
your kit manual).

---

**Q: How many replicates do I need?**

For QC and visual inspection: 1 is sufficient.
For HiCRep concordance: 2+ (score reported as `NOT_ASSESSED` for single replicates).
For differential analysis: ≥ 3 per group (DESeq2 requires at least 2).

---

**Q: What is the difference between `HiCRep NOT_ASSESSED` and `PASS`?**

`NOT_ASSESSED` means there is only one replicate for this sample's subject+mark
combination. We cannot assess whether the library is reproducible because there is
nothing to compare it against. This is **not the same as PASS** — do not interpret
it as validation.

---

**Q: Can I run only the ORACLE export step on an existing `.mcool` and loop BED?**

Yes. Populate `results/cool/<sample>.mcool`, `results/peaks/<sample>_peaks.bed`,
`results/loops/<sample>/*.interactions_*.bed`, and the insulation/eigenvector files,
then run:

```bash
snakemake -s workflow/Snakefile --cores 8 export_oracle --configfile config/config.yaml
```

---

**Q: The `.pt` output file is 0 bytes. What happened?**

PyTorch or PyTorch Geometric failed to import (usually a CUDA library mismatch).
The pipeline now fails the export instead of writing an empty `.pt`.
Check the log: `results/logs/export_oracle/<sample>.log`.

For CPU-only machines, remove `pytorch-cuda` from `environment.yml` and rebuild
the environment with the CPU PyTorch build.

---

**Q: Where do I put the ORACLE export files for the foundation model?**

The `.pt` files are the direct input to the ORACLE model-development corpus.
Their location and provenance are documented in `<sample>.manifest.json`.
Consult the ORACLE model repository (`oracle-model`) for the expected directory
structure and how to register new samples.
