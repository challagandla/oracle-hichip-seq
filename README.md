# ORACLE HiChIP-Seq analysis pipeline

End-to-end HiChIP processing; paired-end FASTQ in, ORACLE-ready Chromatin Operating System (COS) graph prototypes out, with comprehensive QC and publication-grade visualisation at every stage.

> Component of the **ORACLE** research program (Onco-Regulatory Architecture and Chromatin Latent Engine). The proposal and foundation-model training code live in separate repositories.

📖 **New here? Start with the [step-by-step tutorial →](TUTORIAL.md)**

## Why these tool choices

| Choice | What we use | Why |
|---|---|---|
| Workflow manager | **Snakemake 8.x** | Native Python, deterministic DAGs, easy to inspect and extend. |
| Environment | **mamba / micromamba + conda-forge + bioconda + pytorch/nvidia/pyg** | Reproducible `environment.yml`; explicit GPU channels for PyTorch/PyG. |
| Alignment + pairs | **`bwa-mem2 -SP5M` + `pairtools` (Open2C)** | Modern HiC/HiChIP standard, replaces legacy HiC-Pro; emits `.pairs.gz` consumable by `cooler`. |
| Storage | **`.cool` / `.mcool` (cooler)** | Compressed, HDF5-backed, multi-resolution. Generates `.hic` only on demand for Juicebox. |
| Loop calling | **FitHiChIP `Peak-to-ALL`** | Peak-anchored HiChIP loop caller; models 1D ChIP bias. Mustache used as cross-check. |
| Anchor calling | **MACS3** on single read ends | Narrow for K27ac/K4me3/CTCF/cohesin; broad for K27me3/K36me2/K36me3 **and K4me1** (ENCODE). Read ends are ligation partners, not fragment ends, so anchors are called from ends taken individually (`-f BED --nomodel --extsize 147`), never with `-f BAMPE`. MACS2 is not usable: its bioconda binaries link `__log_finite`, dropped from glibc at 2.27. |
| Architectural stripes | **stripenn** at 10 kb | A stripe is a continuous line of enrichment anchored at one end; a loop caller tests discrete pixel pairs and cannot find one. Interpretable as an extrusion anchor on CTCF/cohesin; weaker on H3K27ac, where it is reported and not leaned on. |
| Replicate QC | **HiCRep (`hicrep`)** | Stratum-adjusted correlation robust to distance decay. Grouped on cell type + mark, not donor. SCC is depth-dominated — HiCRep downsamples the deeper library to the shallower — so pairs below `hicrep.min_contacts_for_scc` are flagged `depth_confounded` and cannot decide PASS/FAIL. Missing comparisons report `NOT_ASSESSED`, not PASS. |
| QC suite | **`cooltools` + `pairtools stats` + `MultiQC`** | Cis/trans, distance-decay, P(s), insulation, compartments, APA. |
| Differential loops | **`pyDESeq2`** on per-loop counts | Comparison definitions are explicit and guarded against mixing marks/tissues/protocols. |
| Visualisation | **`pyGenomeTracks` + cooltools APA + HiGlass optional** | Publication-grade static figures; interactive browsing for collaborators. |
| ORACLE export | Custom — multi-resolution graph + prototype node features → `.h5` / PyTorch Geometric `.pt` | Feeds directly into the ORACLE model-development corpus. Current node signal is peak-overlap prototype; true continuous per-mark tracks should be merged later from sister modality pipelines. |

## Layout

```
oracle-hichip/
├── Snakefile                         # root wrapper → workflow/Snakefile
├── README.md, LICENSE, CITATION.cff
├── environment.yml, setup_env.sh
├── config/
│   ├── config.yaml                   # pipeline parameters
│   ├── samples.tsv                   # one row per sample
│   └── genome.yaml                   # genome assemblies + paths
├── workflow/
│   ├── Snakefile                     # main entry point
│   ├── rules/
│   │   ├── 01_qc_raw.smk
│   │   ├── 02_align_pairs.smk
│   │   ├── 03_cool_matrix.smk
│   │   ├── 04_peaks.smk
│   │   ├── 05_loops_fithichip.smk
│   │   ├── 06_loop_qc.smk
│   │   ├── 07_differential.smk
│   │   ├── 08_viz.smk
│   │   ├── 09_export_oracle.smk
│   │   └── 10_multiqc.smk
│   └── scripts/
│       ├── cooltools_eigs_cis.py
│       ├── bedpe_annotate.py
│       ├── export_oracle_cos.py
│       ├── differential_loops.py
│       └── ...
└── docs/
    ├── BEST_PRACTICES.md
    └── ORACLE_INTEGRATION.md
```

## Install

```bash
bash setup_env.sh
mamba activate oracle-hichip
```

For CPU-only machines, remove `pytorch-cuda` from `environment.yml` or maintain a separate CPU environment file.

## Run

```bash
mamba activate oracle-hichip

# Dry run / inspect DAG from repository root
snakemake -s workflow/Snakefile -n --configfile config/config.yaml

# Full run, 32 cores
snakemake -s workflow/Snakefile --cores 32 --configfile config/config.yaml --use-conda

# On a SLURM cluster, provide your own profile
snakemake -s workflow/Snakefile --profile profiles/slurm --configfile config/config.yaml
```

Convenience targets:

```bash
snakemake -s workflow/Snakefile --cores 16 qc_raw --configfile config/config.yaml
snakemake -s workflow/Snakefile --cores 32 align_pairs --configfile config/config.yaml
snakemake -s workflow/Snakefile --cores 32 loops_fithichip --configfile config/config.yaml
snakemake -s workflow/Snakefile --cores 16 export_oracle --configfile config/config.yaml
snakemake -s workflow/Snakefile --cores 4  multiqc --configfile config/config.yaml
```

## Crucial steps and rationale

### 1. QC the raw reads — never skip
FastQC + `fastp` adapter/quality trim. HiChIP libraries often carry adapter readthrough on short fragments; untrimmed reads inflate junk pairs and depress valid-pair yield. Target ≥ Q20 and ≥ 70% retained after trim.

### 2. Align with bwa-mem2 `-SP5M`
`-SP5M` is the canonical HiC/HiChIP mode. Read IDs are now retained in pairtools output so the downstream FitHiChIP validPairs conversion has a valid first column. Only `UU` pairs are retained for loop/matrix work.

### 3. Dedup with `pairtools dedup`, not Picard
HiC duplicates need to be defined on the pair, not either read alone. Expect 10–40% duplicate rate; >50% indicates an undersampled library.

### 4. Build `.cool` + `.mcool` at standard resolutions
The pipeline balances 5 / 10 / 25 / 50 / 100 / 250 / 500 kb / 1 / 2.5 Mb. ORACLE consumes 5 kb / 25 kb / 100 kb / 1 Mb.

### 5. Call 1D peaks on the same reads
HiChIP loop calling needs ChIP anchors. MACS3 mode is chosen per mark in `config/config.yaml`.

### 6. Loop calling — FitHiChIP
FitHiChIP receives a validPairs file generated from pairtools `.pairs.gz` with read IDs preserved. FitHiChIP numeric settings are explicit in `config.yaml`, so changing thresholds does not silently desynchronise the expected output path.

### 7. QC the loops
- Cis/trans ratio ≥ 70% cis.
- Reads in loops ≥ 10% of valid intra-chromosomal pairs.
- P(s) should show expected distance decay.
- APA ≥ 1.5 against a random-shift control. The control moves *both* anchors by one
  offset, preserving the loop's genomic separation, so it is distance-matched by
  construction. Window corners are not: a pixel `(i, j)` sits at separation
  `D + (j - i) * bin`, so corners on opposite sides of the anti-diagonal are not
  comparable to the centre.
- HiCRep SCC as a sanity floor only (`hicrep.threshold_pass`), grouped on cell type
  + mark. SCC is depth-dominated, so a pair whose shallower member falls below
  `hicrep.min_contacts_for_scc` is flagged `depth_confounded` and cannot decide
  PASS/FAIL. Do not read an absolute SCC as evidence that replicates are tighter
  than conditions — on a shallow library it is not.
- Single-replicate HiCRep is reported as `NOT_ASSESSED`, not pass.

### 8. Differential loops
Differential analysis is disabled by default until matched, same-mark case/control groups are defined in `config.yaml`. This prevents accidental comparisons such as tumor H3K27ac versus healthy PBMC H3K36me2.

Example:

```yaml
differential:
  comparisons:
    - name: ovarian_tumor_vs_adjacent_H3K27ac
      mark: H3K27ac
      case_filter:    { tissue: tumor, disease: ovarian_HGSOC, library_protocol: HiChIP_v2 }
      control_filter: { tissue: adjacent_normal, disease: non_cancer, library_protocol: HiChIP_v2 }
```

### 9. Architectural stripes
Called with stripenn at 10 kb, not at the 5 kb loop resolution: a stripe is detected
as an image feature (Canny edge detection over the contact map), and at 5 kb a
HiChIP matrix is sparse enough that the edges traced are mostly sampling noise.

Stripes are not recoverable from a loop list. A loop caller tests discrete pixel
pairs against a distance-decay background, so a stripe — a continuous line of
enrichment running away from one anchor — is at best fragmented into a row of
individually unconvincing pixels. Read them against the anchor: on CTCF or cohesin a
stripe is directly a loop-extrusion anchor; on H3K27ac the anchors are enhancers,
extrusion is not what defines them, and fewer stripes is the expected result rather
than a failure. Reports split by mark and never pool them.

### 10. Visualisation
- `pyGenomeTracks` for arc + heatmap composite figures.
- APA for aggregate loop strength.
- Virtual 4C from anchor of interest.
- Five cohort figures in `results/figures/` (vector PDF + 400 dpi PNG), including
  loop yield plotted *against* library depth — the panel that distinguishes a
  biological difference in loop count from a sequencing one. Libraries below the
  depth floor are marked, not dropped.

### 11. Export to ORACLE COS format
Per sample the exporter emits:

- `cos_<sample>.h5` equivalent at `results/oracle_cos/<sample>.h5`
- PyTorch Geometric graph at `results/oracle_cos/<sample>.pt`
- manifest JSON describing feature channels and limitations

Current node features are `peak_overlap_count_per_kb`, insulation and E1 eigenvector. Treat this as a structural/peak prototype, not the final full multimodal COS.

## License

MIT — see [`LICENSE`](LICENSE).

## Contact

Anil Challagandla — challagandla.anil@gmail.com
