# ORACLE prototype export contract

Stage 09 converts each HiChIP library into a multi-resolution structural graph for
downstream ORACLE development. This is a prototype, HiChIP-only representation. It
is not yet a complete multimodal model input and should not be described as a
finished or canonical training schema.

## Files

For each sample, the workflow writes:

- `results/oracle_cos/<sample>.pt` — PyTorch Geometric `HeteroData`;
- `results/oracle_cos/<sample>.h5` — HDF5 representation with genomic coordinates;
- `results/oracle_cos/<sample>.manifest.json` — schema and sample metadata, channel
  order, feature availability and normalization, resolutions, blacklist accounting,
  declared-input hashes, output hashes, and scope note.

The rule emits both `.pt` and `.h5` so tensor consumers and language-independent
inspection use the same graph content.

Only load `.pt` files produced by a trusted workflow run. PyTorch serialization is
not a safe exchange format for untrusted files. The HDF5 file is the better surface
for language-independent inspection.

## Resolutions and node order

The default graph resolutions are 5 kb, 25 kb, 100 kb, and 1 Mb, configured by
`cooler.oracle_bin_sizes_kb` in `config/config.yaml`.

The schema identifier is `oracle-hichip-cos-v2`, and all coordinates are 0-based,
half-open. Each HDF5 group is named `res_<bin_size_bp>` and contains:

- `bin_chrom`, `bin_start`, and `bin_end`;
- `x`, `x_observed_mask`, and `blacklist_mask` for node features and validity;
- `loop_edge_index` and `loop_edge_attr`;
- `adjacency_edge_index`;
- `n_nodes`, `n_loop_edges`, `n_adjacency_edges`, and `bin_size_bp` attributes.

The HDF5 root records `node_feature_channels` and
`loop_edge_attr_channels`; the manifest and PyG object use the same names and
ordering.

The HDF5 `hierarchy/res_<fine>_to_res_<coarse>` groups contain both
`contained_by_edge_index` and its reversed `contains_edge_index`. Successive ORACLE
resolutions must divide exactly, and every fine node has exactly one parent at the
next coarser resolution.

The corresponding PyG node type is `bin_res_<bin_size_bp>`. Each node store carries
`chrom_id`, the ordered `chrom_names` lookup, `bin_start`, `bin_end`,
`x_observed_mask`, and `blacklist_mask`, in addition to `x`. The PyG sample store
carries one explicit sample node (`num_nodes=1`) with sample ID, assembly, mark,
cell type, and optional token tensors. HDF5
root attributes and the manifest record the same sample/schema metadata.
By default `oracle_export.primary_chromosomes_only: true` retains chr1-22/X and
excludes alternate, decoy, unplaced, chrY, and mitochondrial sequences so graph
nodes match the workflow's primary-chromosome structural analyses. If that policy is
disabled, chromosomes listed in `oracle_export.drop_chromosomes` are omitted; by
default these are chrM and chrY.

## Node features

Columns in `x`, in order, are:

1. `peak_overlap_count_per_kb` — number of overlapping per-sample MACS3 peak
   intervals, normalized by bin length in kilobases;
2. `insulation` — the configured local insulation track projected by interval
   overlap;
3. `E1_eigenvector` — the GC-phased exploratory E1 track projected by interval
   overlap.

Coarse tracks are broadcast to every finer bin they cover, while coarser destination
bins receive an overlap-weighted mean. Bins without insulation or E1 coverage are
filled with zero to preserve tensor shape, and `x_observed_mask` is false for those
channel/bin entries. `blacklist_mask` separately identifies affected nodes.

The channel order does not change when a feature is unavailable. Instead, the
manifest contains a `node_feature_availability` object keyed by the same channel
names. Each entry records `available` and `normalization`, plus balance status and
resolution where balance is relevant. In particular:

- `insulation` remains available after either converged ICE or the documented
  raw-count fallback; its normalization tells consumers which one was used;
- `E1_eigenvector` is `available: false` when balancing did not pass at 100 kb, and
  its entry records why;
- zeros retained in an unavailable channel preserve tensor shape only. They are
  missing-value fill, not observed biological zeros.

Consumers must combine the fixed column order with `node_feature_availability`,
`x_observed_mask`, and `blacklist_mask`; inspecting `x` alone cannot distinguish an
unavailable feature from a true numeric zero.

Bins overlapping the configured blacklist have every feature cell marked
unobserved and zero-filled, not merely the peak channel. The manifest records
whether blacklist masking was applied, the blacklist hash, and the number of
affected bins at each resolution.

`peak_overlap_count_per_kb` is an interval-overlap prototype, not continuous
per-mark signal. A full multimodal representation should add measured and
normalized signal channels from appropriate ChIP-seq, CUT&Tag, ATAC-seq,
methylation, and expression assays.

## Edges

Every resolution contains two separate directed relations:

- `adjacent`: consecutive genomic bins on the same chromosome, in both directions;
- `loop`: filtered FitHiChIP contacts, in both directions.

Only `loop` has attributes. `loop_edge_attr` columns, in order, are:

1. `max_loop_score` — maximum FitHiChIP contact-count/score among fine calls that
   map to the node pair;
2. `min_loop_pvalue` — minimum raw P value among those calls;
3. `min_loop_fdr` — minimum adjusted Q/FDR value among those calls;
4. `genomic_distance_bp`;
5. `fine_loop_count` — number of fine calls consolidated into this coarse pair.

Parallel fine calls are consolidated so generic message passing does not weight a
coarse pair by accidental row multiplicity. Minimum P/Q are strongest-evidence
summaries, not newly tested coarse-bin significance values. Contacts whose anchors
collapse into the same coarse node are omitted. Nearest-gene and peak annotations
used in the intermediate BEDPE are not currently exported as edge attributes.

Successive node types are connected by `contained_by` (fine to coarse) and
`contains` (coarse to fine) relations. This makes the four scales one explicit
hierarchy rather than four disconnected graphs. The manifest enumerates the exact
PyG edge-type triplets for every within-resolution and hierarchy relation.

## Optional sample-level tokens

Set `oracle_export.microbiome_metadata_tsv` explicitly to enable sample-level
tokens; the default is empty. A configured table is a declared Snakemake input, so
creating or changing it invalidates the export, and its SHA-256 is added to the
manifest. The table must contain a unique `sample_id` column. Columns containing
numeric observations define one table-wide key order. Values and their observed
mask are stored with shape `[1, F]` on the one sample node and mirrored as
`microbiome_values` and `microbiome_observed_mask` in HDF5. A missing cell or absent
sample is zero-filled with a false mask, so all samples retain the same batchable
width; a column that mixes numeric and non-numeric observed values fails rather than
silently changing schema. Entirely non-numeric annotation columns are skipped. A
configured file that is missing fails dependency resolution.

## QC and provenance

The export depends on the per-sample QC, balance, and E1-status evidence at:

```text
results/qc/loop_qc/<sample>.json
results/qc/balance/<sample>.balance.json
results/qc/compartments/<sample>.cis.eigs.status.json
```

Balance is assessed per resolution:

- `PASS`: the ICE weight exists and cooler reports that it converged;
- `WARN`: the balancing attempt did not converge, so the attempted weight is not
  published or used;
- `NOT_ASSESSED`: the weight or convergence evidence is absent.

`WARN` does not invalidate raw contacts, and `NOT_ASSESSED` is not a pass. The
manifest carries the balance summary and its hash so an ORACLE consumer can audit
the normalization decision without reopening the cooler.

Other component states include:

- contact-map-held-out APA: `DESCRIPTIVE` or `NOT_ASSESSED`;
- HiCRep: `PASS`, `FAIL`, `DISCORDANT`, or `NOT_ASSESSED`;
- overall QC: `PASS`, `FAIL`, `PASS_WITH_UNCERTAINTY`, or
  `PASS_WITH_NOT_ASSESSED`.

`NOT_ASSESSED` is not a pass. Downstream ingestion should decide explicitly whether
to exclude, retain, or mask samples with incomplete evidence.
`DISCORDANT` maps to `PASS_WITH_UNCERTAINTY`, never a machine pass
(`overall_pass: false`), unless another hard QC failure makes the overall state
`FAIL`.

The per-sample manifest records schema, assembly, mark, cell type, coordinate system,
node-order description, `node_feature_availability`, the balance summary, and hashes
for the cooler, balance report, annotated loop table, peaks, insulation, E1 and its
status, loop-QC input, chromosome sizes, blacklist, export helper code, `.pt`
output, and HDF5 output. The separate
`results/provenance/run_manifest.json` embeds and hashes the effective pipeline
configuration, sample sheet, genome configuration, actual reference files, resolved
packages for portable rule environments, QC summaries, differential universe/design
records, report/figure files, and ORACLE manifests, and records the Git revision.
Keep both manifests with an exported dataset.

Subject ID and the full experimental design remain in `config/samples.tsv`, not the
per-sample graph. Preserve the sample sheet, pipeline/genome configuration, and
HDF5 companion as part of the dataset contract.

## Differential results are separate

Differential contact tables are not embedded in the per-sample graphs. Comparisons
are configured explicitly, use a condition-blind exact-pixel universe derived from
unthresholded FitHiChIP all-interaction tables, and require an explicit donor/pairing
factor in this release. Read each comparison's
`hypothesis_universe.json`, `design.json`, and `paired_effects.tsv` before using its
results. `PILOT_UNDERPOWERED` output is exploratory even when a Wald adjusted p-value
crosses a conventional threshold. A differential HiChIP signal can reflect anchor
occupancy, contact frequency, or both.

## Minimal downstream checks

Before ingesting an export:

1. verify every manifest hash and the expected Git revision;
2. confirm schema, assembly, sample identity, mark, and cell type against the
   manifest, and subject/design metadata against the retained sample sheet;
3. confirm resolution and feature/edge column order against the manifest;
4. inspect every `node_feature_availability` entry and its normalization before
   reading that column from `x`;
5. inspect balance and component QC states rather than only `overall_pass`;
6. keep missing-value fill, blacklist masking, and dropped chromosomes explicit;
7. validate node/edge counts and coordinate order against the HDF5 file;
8. do not treat nearest-gene annotation as causal enhancer-target evidence.
