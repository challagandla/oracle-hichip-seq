# ORACLE integration contract

This repository exports HiChIP-derived graph objects for the ORACLE model-development corpus.

## Outputs

For each sample:

- `results/oracle_cos/<sample>.pt` — PyTorch Geometric `HeteroData` graph.
- `results/oracle_cos/<sample>.h5` — HDF5 mirror for inspection and non-PyG consumers.
- `results/oracle_cos/<sample>.manifest.json` — machine-readable description of resolutions, node features, edge features and caveats.

## Resolutions

The graph is emitted at:

- 5 kb
- 25 kb
- 100 kb
- 1 Mb

These are configured by `cooler.oracle_bin_sizes_kb` in `config/config.yaml`.

## Current node features

The current HiChIP-only export contains:

1. `peak_overlap_count_per_kb` — MACS3 peak-overlap count normalised by bin length.
2. `insulation` — cooltools insulation score.
3. `E1_eigenvector` — A/B compartment eigenvector.

Important: `peak_overlap_count_per_kb` is a prototype feature, not continuous per-mark ChIP/CUT&Tag signal. Full ORACLE COS should merge continuous signal channels from sister ATAC/CUT&Tag/methylation/RNA pipelines.

## Current edge features

Edges include:

- genomic adjacency edges
- FitHiChIP loop/contact edges

Loop edge attributes:

1. `loop_score`
2. `loop_fdr`
3. `genomic_distance_bp`

## QC contract

Each sample has a QC JSON at:

`results/qc/loop_qc/<sample>.json`

Replicate QC is three-state:

- `PASS`
- `FAIL`
- `NOT_ASSESSED`

Single-replicate samples should not be interpreted as replicate-validated.

## Differential analysis contract

Differential loops must be configured explicitly and stratified by mark/tissue/protocol. The pipeline intentionally leaves `differential.comparisons` empty in the template to avoid invalid default comparisons.
