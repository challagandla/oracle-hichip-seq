# How this pipeline connects to the ORACLE foundation model

The HiChIP pipeline emits one `cos_<sample>.pt` file per sample (PyTorch
Geometric `HeteroData`) plus a parallel `.h5` mirror. These are the canonical
inputs to ORACLE's two-tower training code.

## File contract (output of `09_export_oracle`)

```
results/oracle_cos/<sample>.pt        # PyG HeteroData (primary)
results/oracle_cos/<sample>.h5        # HDF5 mirror (debugging / non-PyG consumers)
results/oracle_cos/<sample>.manifest.json
results/oracle_cos/intermediates/<sample>.annotated_loops.bedpe
```

### HeteroData layout

```
data['sample'].id                              # str
data['sample'].microbiome                      # FloatTensor[K] (may be empty)

for bp in [5000, 25000, 100000, 1_000_000]:
    nt = f'bin_res_{bp}'
    data[nt].x                                 # FloatTensor[N_bp, C]
                                               #   C currently = [peak_signal,
                                               #     insulation, E1] = 3
                                               #   sister pipelines (atac/, cutandtag/,
                                               #     medip/, rnaseq/) append channels
                                               #     during merge in oracle/training/data.py.
    data[nt].bin_chrom / bin_start / bin_end   # stored in .h5; restored on load if needed

    et = (nt, 'contact', nt)
    data[et].edge_index                        # LongTensor[2, E]
    data[et].edge_attr                         # FloatTensor[E, 3]
                                               #   [score, fdr, distance_bp]
    data[et].edge_kind                         # Int8Tensor[E]
                                               #   0 = adjacency, 1 = HiChIP loop
```

### Where additional modality channels come in

| Channel | Produced by | Merged where |
|---|---|---|
| `peak_signal` (HiChIP anchor peaks) | this pipeline | export_oracle_cos.py |
| H3K4me3 / K27me3 / K36me2 / K36me3 CUT&Tag bigWig sums per bin | `code/cutandtag/` | `oracle/training/data.py` |
| ATAC accessibility | `code/atac/` | `oracle/training/data.py` |
| MeDIP CpG methylation | `code/medip/` | `oracle/training/data.py` |
| RNA expression (bin–gene aware) | `code/rnaseq/` | `oracle/training/data.py` |
| WGS SVs → edge rewiring | `code/wgs_sv/` | `oracle/training/data.py` (differentiable rewire layer) |
| Microbiome tokens | `config/microbiome_tokens.tsv` here | this pipeline writes `data['sample'].microbiome` directly |

## Training-side loader (illustrative — lives in `oracle/training/data.py`)

```python
import torch
from torch_geometric.data import HeteroData

def load_sample(path: str) -> HeteroData:
    return torch.load(path, weights_only=False)

# Merge in other modalities (called once per epoch or pre-baked):
def attach_external_channels(data: HeteroData, sample_id: str,
                             cutandtag_h5: str, atac_h5: str,
                             medip_h5: str, rna_h5: str,
                             wgs_sv_h5: str) -> HeteroData:
    ...
    # Concatenate new channels onto data[node_type].x and append
    # SV-induced edges to data[(...,'contact',...)].edge_index.
    return data
```

## Cross-pipeline conventions

- **Bin coordinates are inclusive-exclusive [start, end)**, always at the
  resolutions 5 kb / 25 kb / 100 kb / 1 Mb. Never bin to other sizes downstream.
- **Genome assembly is recorded in the manifest.** Mixing hg38 and T2T-CHM13
  channels in the same `HeteroData` is a hard error.
- **Blacklisted bins are masked (not removed).** Mask flag is stored in the
  sister pipelines' channel exports as an extra "is_blacklist" feature.
- **Channels are z-scored per-modality within a batch** during training
  (not here) — this pipeline emits raw per-bin signal.

## QC gate before a sample enters training

A sample's `.pt` is only consumed by training if its
`results/qc/loop_qc/<sample>.json` carries `overall_pass: true`. The
ORACLE training script (`oracle/training/build_dataset.py`) reads the QC
JSONs and discards failing samples with a logged reason.
