# ORACLE HiChIP analysis pipeline

End-to-end HiChIP processing — paired-end FASTQ in, ORACLE-ready Chromatin Operating System (COS) graphs out, with comprehensive QC and publication-grade visualisation at every stage.

> Component of the **ORACLE** research program (Onco-Regulatory Architecture and Chromatin Latent Engine). The proposal and foundation-model training code live in separate repositories.

## Why these tool choices

| Choice | What we use | Why |
|---|---|---|
| Workflow manager | **Snakemake 8.x** | Native Python, deterministic DAGs, trivial to step through in VS Code, easy to extend per-sample. |
| Environment | **mamba / micromamba + conda-forge + bioconda** | Solver 10–20× faster than vanilla conda; reproducible `environment.yml`. |
| Alignment + pairs | **`bwa-mem2 -SP5M` + `pairtools` (Open2C)** | Modern standard, replaces legacy HiC-Pro; outputs `.pairs.gz` consumable by `cooler`. |
| Storage | **`.cool` / `.mcool` (cooler)** | Compressed, HDF5-backed, multi-resolution. Generates `.hic` only on demand for Juicebox. |
| Loop calling | **FitHiChIP `Peak-to-ALL`** | Gold standard for peak-anchored HiChIP; correctly models the 1D ChIP bias. Mustache used as cross-check. |
| Peak calling | **MACS2** (narrow for K27ac/K4me3, broad for K27me3/K36me2/me3) | Standard; outputs feed FitHiChIP anchors. |
| Replicate QC | **HiCRep (`hicreppy`)** | Stratum-adjusted correlation — only metric robust to distance decay. |
| QC suite | **`cooltools` + `pairtools stats` + `MultiQC`** | Cis/trans, distance-decay, P(s), insulation. |
| Normalisation | **ICE (cooler balance) + KR** | ICE for matrices; loops use FitHiChIP spline correction internally. |
| Differential loops | **`pyDESeq2`** on per-loop counts (alt: `diffHic` R) | Mature framework with optional R fallback. |
| Visualisation | **`pyGenomeTracks` + cooltools APA + HiGlass (optional)** | Publication-grade static; interactive browsing for collaborators. |
| ORACLE export | Custom — multi-resolution graph + node features → `.h5` / PyTorch Geometric `.pt` | Feeds directly into the ORACLE foundation-model training corpus. |
| IDE | **VS Code** with Snakemake, Python, Jupyter, Pylance extensions | See `.vscode/extensions.json`. |

## Layout

```
oracle-hichip/
├── README.md                       # this file
├── LICENSE                         # MIT
├── CITATION.cff
├── .gitignore, .gitattributes
├── environment.yml                 # conda/mamba env
├── setup_env.sh                    # one-command install
├── config/
│   ├── config.yaml                 # pipeline parameters
│   ├── samples.tsv                 # sample sheet (one row per sample)
│   └── genome.yaml                 # genome assemblies + paths
├── workflow/
│   ├── Snakefile                   # main entry point
│   ├── rules/
│   │   ├── 01_qc_raw.smk           # FastQC + fastp adapter trim
│   │   ├── 02_align_pairs.smk      # bwa-mem2 + pairtools parse/sort/dedup
│   │   ├── 03_cool_matrix.smk      # cooler cload + zoomify + balance (ICE)
│   │   ├── 04_peaks.smk            # MACS2 from 1D reads
│   │   ├── 05_loops_fithichip.smk  # pairs→validpairs + FitHiChIP + mustache
│   │   ├── 06_loop_qc.smk          # cis/trans, P(s), insulation, eigs, APA, HiCRep
│   │   ├── 07_differential.smk     # union loops + pyDESeq2 differential
│   │   ├── 08_viz.smk              # pyGenomeTracks + virtual 4C
│   │   ├── 09_export_oracle.smk    # COS graph + node features → ORACLE
│   │   └── 10_multiqc.smk          # aggregate HTML report
│   └── scripts/                    # Python: HiCRep, APA, loop QC, BEDPE annotate,
│                                   # differential loops, pyGenomeTracks, virtual 4C,
│                                   # ORACLE COS exporter, utils, etc.
├── .vscode/                        # VS Code settings + recommended extensions
└── docs/
    ├── BEST_PRACTICES.md           # 46-item analysis checklist
    └── ORACLE_INTEGRATION.md       # contract with the ORACLE foundation model
```

## Install

```bash
# install micromamba (fastest) — skip if you have mamba/conda already
curl -L micro.mamba.pm/install.sh | bash

# build the env
mamba env create -f environment.yml -n oracle-hichip
mamba activate oracle-hichip

# verify
snakemake --version && cooler --version && pairtools --version && macs2 --version
```

## Run

```bash
mamba activate oracle-hichip

# dry run to inspect the DAG
snakemake -n --configfile config/config.yaml

# full run, 32 cores
snakemake --cores 32 --configfile config/config.yaml --use-conda

# on a SLURM cluster (provide your own profile)
snakemake --profile profiles/slurm --configfile config/config.yaml
```

Convenience phony targets:

```bash
snakemake --cores 16 qc_raw
snakemake --cores 32 align_pairs
snakemake --cores 32 loops_fithichip
snakemake --cores 16 export_oracle
snakemake --cores 4  multiqc
```

## Crucial steps and rationale

### 1. QC the raw reads — never skip
FastQC + `fastp` adapter/quality trim. HiChIP libraries often carry adapter readthrough on short fragments; untrimmed reads inflate junk pairs and depress valid-pair yield. Target ≥ Q20 and ≥ 70% retained after trim.

### 2. Align with bwa-mem2 `-SP5M`
`-SP5M` is the canonical HiC/HiChIP mode (skip mate rescue, soft-clip 5′ supplementary, mark short hits as secondary). pairtools parses the alignment into a `.pairs.gz` file with the four-letter pair type code. **Only `UU` pairs are uniquely-mapped and trustworthy.**

### 3. Dedup with `pairtools dedup` (NOT Picard)
HiC duplicates need to be defined on the pair (read1 position, read2 position) not on either read alone. Picard MarkDuplicates undercounts or overcounts. Expect 10–40% duplicate rate; > 50% indicates undersampled library.

### 4. Build `.cool` + `.mcool` at the standard resolutions
We balance at 5 / 10 / 25 / 50 / 100 / 250 / 500 kb / 1 / 2.5 Mb. ORACLE consumes 5 kb / 25 kb / 100 kb / 1 Mb as the four hierarchy levels.

### 5. Call 1D peaks on the same reads (MACS2)
HiChIP loop calling needs ChIP anchors. We extract single-end reads from the deduped UU pairsam (`pairtools split`) and run MACS2:
- H3K27ac / H3K4me3 / H3K4me1 / CTCF → **narrowPeak** (q < 0.01)
- H3K27me3 / H3K36me2 / H3K36me3 → **broadPeak** (q < 0.05)

### 6. Loop calling — FitHiChIP (Peak-to-ALL)
At 5 kb bin size, FDR < 0.01, ≥ 6 reads per loop, 20 kb–3 Mb distance range. Always report loops at each resolution separately — collapsing across resolutions inflates FDR.

### 7. QC the loops
- **Cis/trans ratio** ≥ 70% cis.
- **Reads in loops** ≥ 10% of valid intra-chromosomal pairs.
- **Distance decay P(s)** should follow the expected −1 slope on log–log.
- **APA score** ≥ 1.5 at high-confidence loops vs. random shifts.
- **HiCRep stratum-adjusted correlation** between biological replicates ≥ 0.85.

### 8. Differential loops
Union loop set across samples, count per-loop per sample with `cooler.matrix(...).fetch`, run `pyDESeq2` with FDR < 0.05 and |log2FC| ≥ 1.

### 9. Visualisation
- `pyGenomeTracks` for arc + heatmap composite figures.
- Aggregate Peak Analysis (APA) for loop strength.
- Virtual 4C from anchor of interest.
- Optional: HiGlass server for interactive sharing.

### 10. Export to ORACLE COS format
Per sample we emit:
- `cos_<sample>.h5` — node features (signal per bin) + edge list (loops, adjacency) at four resolutions.
- `cos_<sample>.pt` — PyTorch Geometric `HeteroData` covering 5 kb / 25 kb / 100 kb / 1 Mb resolutions.

The exporter in `workflow/scripts/export_oracle_cos.py` is the contract between this pipeline and the ORACLE foundation model. See `docs/ORACLE_INTEGRATION.md`.

## Best practices

See [`docs/BEST_PRACTICES.md`](docs/BEST_PRACTICES.md) — a 46-item checklist covering library design, alignment, dedup, matrix balance, peak/loop calling, replicate concordance, differential analysis, visualisation, and reproducibility.

## License

MIT — see [`LICENSE`](LICENSE).

## Contact

Anil Challagandla — challagandla.anil@gmail.com
