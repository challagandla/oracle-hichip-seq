# Best practices for ORACLE HiChIP analysis

A condensed checklist. Treat each numbered item as a non-negotiable unless explicitly justified in the project log.

## Library design

1. **Always sequence paired-end ≥ 75 bp.** Shorter reads collapse cis/trans calling.
2. **Three biological replicates per condition is the minimum** for differential analysis. Do not accept them on an absolute HiCRep SCC: HiCRep downsamples the deeper library to the shallower, so SCC grades the worst member's depth. Measured on GSE101498 at 25 kb, two libraries from *different* cell types (90.8M and 47.6M contacts) score 0.743 while two genuine replicates score 0.221 because one holds 3.0M — an absolute cutoff fails the real replicate pair and passes the cross-cell-type one. Sequence every replicate past `hicrep.min_contacts_for_scc` and only then read the SCC.
3. **Depth is the first thing to rule out, not the last.** Loop count, stripe count and SCC all rise with unique cis pairs. Before reporting any of them as a difference between conditions, plot it against library depth (`results/figures/figure3_loops_apa`).
4. **Match library prep batch across conditions** when possible. Otherwise treat batch as a covariate in differential analysis.

## Read QC

5. Always run FastQC + fastp. Discard libraries with retention < 70% after Q20 + length 36 filter.
6. **Never analyse reads that haven't passed adapter trimming.** HiChIP libraries often have adapter readthrough; pairtools will silently miscall pair types on adapter-contaminated reads.

## Alignment

7. Use `bwa-mem -SP5M` (or `bwa-mem2 -SP5M`). Never use bowtie2 or BWA-aln for HiChIP.
8. Use pairtools, not HiC-Pro. pairtools is actively maintained, faster, and produces the canonical `.pairs.gz` format that every downstream tool consumes.
9. **Keep only `UU` (unique-unique) pair types.** Anything else introduces ambiguity and inflates false-positive loops.
10. Set `--min-mapq 30`. Lower values let multi-mappers contaminate the contact map.

## Deduplication

11. **Always dedup with pairtools, never with Picard MarkDuplicates.** Duplicate detection must use both ends of the pair.
12. Library complexity check: > 50% duplicate rate → library is undersampled. Either resequence or flag.
13. Valid pair yield (UU after dedup / total input pairs) ≥ 25% is the floor for further analysis.

## Matrix

14. Build `.mcool` at standard resolutions (5/10/25/50/100/250/500 kb, 1/2.5 Mb). ICE-balance every resolution.
15. **Never analyse balanced and unbalanced matrices interchangeably.** Loop calling uses raw counts; matrix similarity uses balanced.
16. Use cooler/cool format, not legacy `.hic`. Convert to `.hic` only for Juicebox display.

## Distance decay (P(s))

17. The expected log–log P(s) slope is approximately −1 (in cis). Deviations:
    - Slope shallower than −0.8 → unusual long-range bias (over-amplified or chimeric library).
    - Slope steeper than −1.3 → fragment-length bias; consider re-sonicating.
18. **Cis/trans ratio ≥ 70%** is the QC floor; < 60% means trans contamination or undersampled.

## 1D peak calling

19. Call peaks from `pairtools split` 1D reads, not from a separate ChIP-seq aliquot — the protein-of-interest immunoprecipitation is implicit in HiChIP.
20. Use narrow mode for punctate marks (H3K27ac, H3K4me3, CTCF, cohesin), broad mode for spreading marks (H3K27me3, H3K36me2, H3K36me3) — and for **H3K4me1**, which is broad per ENCODE despite marking enhancers: its signal is a wide bimodal shoulder around the nucleosome-depleted region, not a sharp peak. Getting this wrong matters more in HiChIP than in ChIP-seq: anchors are what loops are hung on, so calling a punctate mark with `--broad` fuses adjacent enhancers into multi-kb blocks and merges distinct enhancer–promoter contacts into one uninterpretable loop.
21. Q-value 0.01 for narrow, 0.05 for broad. Document any deviation.

## Loop calling

22. Use **FitHiChIP `Peak-to-ALL`** as the primary caller. It correctly models the 1D ChIP bias.
23. Distance range: **20 kb to 3 Mb**. Loops shorter than 20 kb confound with non-loop short-range contacts; > 3 Mb is mostly noise.
24. **Run loop calling per resolution** (5 kb and 10 kb) and report each separately. Never merge across resolutions before FDR control.
25. ≥ 6 supporting reads per loop. Less is anecdotal.
26. FDR threshold q < 0.01.
27. Cross-check with `mustache` at 10 kb. Concordant calls are more trustworthy.
28. **Minimum acceptable loop count is ~1,000 per sample** for FitHiChIP at q < 0.01. Fewer indicates poor library or wrong protein.

## APA

29. APA score ≥ 1.5 vs random-shift controls is the QC pass threshold.
30. Always use distance-stratified random controls — uniform random controls underestimate background.

## Replicate concordance

31. **Group replicates on condition + protein, never on donor.** With several conditions taken from the same donors, grouping on donor puts different cell types into one "replicate" group and reports the biology you are trying to detect as if it were noise.
32. **Do not accept or reject a replicate on an absolute HiCRep SCC.** It is depth-dominated (see item 2). Use it as a sanity floor that the matrix looks like a contact map, check the pair cleared `hicrep.min_contacts_for_scc`, and treat any pair below that floor as `depth_confounded` — reported, but not evidence either way.
33. A low SCC on an adequately-sequenced pair is a real failure: re-prep. A low SCC on a shallow pair is a sequencing decision, not a biology one.

## Differential analysis

34. Build a union loop set across all samples, then count per-loop contacts per sample.
35. Use DESeq2 (or pyDESeq2) with `library size` from valid-pair count, not from raw read count.
36. FDR < 0.05 and |log2FC| ≥ 1 is the standard cutoff.
37. **Always run a permutation test** (shuffle group labels, recompute) to estimate empirical FDR — DESeq2's assumption violations can be subtle in HiChIP count data.

## Visualisation

38. Use pyGenomeTracks for static figures. Always include axis ticks, scale bars, gene model, and resolution annotation.
39. APA plots: use log1p colour scale (linear scale washes out the centre).
40. Virtual 4C: only valid if the viewpoint anchor has > 50 supporting reads at the chosen resolution.

## Reproducibility

41. **Lock the conda environment** after install: `mamba env export --no-builds > environment.lock.yml`.
42. **Pin a random seed** in every script that subsamples (HiCRep, APA controls, differential permutations).
43. Use Snakemake `--use-conda` (or `--use-singularity`) so each rule runs in an isolated env.
44. Tag all intermediate files with the assembly version (`hg38_`, `t2t_`, etc.) in the filename.
45. Commit `samples.tsv`, `config/*.yaml`, `environment.yml`, `environment.lock.yml`, and the Snakefile to git. Do not commit results.

## ORACLE integration

46. Always run `09_export_oracle` after loop calling completes. The `.pt` file is the contract with the foundation model.
47. The ORACLE export carries provenance: the `manifest.json` records every input file's SHA-256 hash. Never modify a `.pt` by hand.
48. Sequence-imputed marks must be flagged in the `manifest.json`; the training code can then mask them differently from measured marks.
