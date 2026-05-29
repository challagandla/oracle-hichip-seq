# Best practices for ORACLE HiChIP analysis

A condensed checklist. Treat each numbered item as a non-negotiable unless explicitly justified in the project log.

## Library design

1. **Always sequence paired-end ≥ 75 bp.** Shorter reads collapse cis/trans calling.
2. **Three biological replicates per condition is the minimum** for differential analysis. HiCRep SCC ≥ 0.85 between any two replicates of the same condition.
3. **Match library prep batch across conditions** when possible. Otherwise treat batch as a covariate in differential analysis.

## Read QC

4. Always run FastQC + fastp. Discard libraries with retention < 70% after Q20 + length 36 filter.
5. **Never analyse reads that haven't passed adapter trimming.** HiChIP libraries often have adapter readthrough; pairtools will silently miscall pair types on adapter-contaminated reads.

## Alignment

6. Use `bwa-mem -SP5M` (or `bwa-mem2 -SP5M`). Never use bowtie2 or BWA-aln for HiChIP.
7. Use pairtools, not HiC-Pro. pairtools is actively maintained, faster, and produces the canonical `.pairs.gz` format that every downstream tool consumes.
8. **Keep only `UU` (unique-unique) pair types.** Anything else introduces ambiguity and inflates false-positive loops.
9. Set `--min-mapq 30`. Lower values let multi-mappers contaminate the contact map.

## Deduplication

10. **Always dedup with pairtools, never with Picard MarkDuplicates.** Duplicate detection must use both ends of the pair.
11. Library complexity check: > 50% duplicate rate → library is undersampled. Either resequence or flag.
12. Valid pair yield (UU after dedup / total input pairs) ≥ 25% is the floor for further analysis.

## Matrix

13. Build `.mcool` at standard resolutions (5/10/25/50/100/250/500 kb, 1/2.5 Mb). ICE-balance every resolution.
14. **Never analyse balanced and unbalanced matrices interchangeably.** Loop calling uses raw counts; matrix similarity uses balanced.
15. Use cooler/cool format, not legacy `.hic`. Convert to `.hic` only for Juicebox display.

## Distance decay (P(s))

16. The expected log–log P(s) slope is approximately −1 (in cis). Deviations:
    - Slope shallower than −0.8 → unusual long-range bias (over-amplified or chimeric library).
    - Slope steeper than −1.3 → fragment-length bias; consider re-sonicating.
17. **Cis/trans ratio ≥ 70%** is the QC floor; < 60% means trans contamination or undersampled.

## 1D peak calling

18. Call peaks from `pairtools split` 1D reads, not from a separate ChIP-seq aliquot — the protein-of-interest immunoprecipitation is implicit in HiChIP.
19. Use narrow mode for promoter/enhancer marks (H3K27ac, H3K4me3, H3K4me1, CTCF), broad mode for spreading marks (H3K27me3, H3K36me2, H3K36me3).
20. Q-value 0.01 for narrow, 0.05 for broad. Document any deviation.

## Loop calling

21. Use **FitHiChIP `Peak-to-ALL`** as the primary caller. It correctly models the 1D ChIP bias.
22. Distance range: **20 kb to 3 Mb**. Loops shorter than 20 kb confound with non-loop short-range contacts; > 3 Mb is mostly noise.
23. **Run loop calling per resolution** (5 kb and 10 kb) and report each separately. Never merge across resolutions before FDR control.
24. ≥ 6 supporting reads per loop. Less is anecdotal.
25. FDR threshold q < 0.01.
26. Cross-check with `mustache` at 10 kb. Concordant calls are more trustworthy.
27. **Minimum acceptable loop count is ~1,000 per sample** for FitHiChIP at q < 0.01. Fewer indicates poor library or wrong protein.

## APA

28. APA score ≥ 1.5 vs random-shift controls is the QC pass threshold.
29. Always use distance-stratified random controls — uniform random controls underestimate background.

## Replicate concordance

30. **HiCRep stratum-adjusted correlation ≥ 0.85** between biological replicates of the same condition and protein.
31. If SCC < 0.7 between replicates, treat one as failed and re-prep.

## Differential analysis

32. Build a union loop set across all samples, then count per-loop contacts per sample.
33. Use DESeq2 (or pyDESeq2) with `library size` from valid-pair count, not from raw read count.
34. FDR < 0.05 and |log2FC| ≥ 1 is the standard cutoff.
35. **Always run a permutation test** (shuffle group labels, recompute) to estimate empirical FDR — DESeq2's assumption violations can be subtle in HiChIP count data.

## Visualisation

36. Use pyGenomeTracks for static figures. Always include axis ticks, scale bars, gene model, and resolution annotation.
37. APA plots: use log1p colour scale (linear scale washes out the centre).
38. Virtual 4C: only valid if the viewpoint anchor has > 50 supporting reads at the chosen resolution.

## Reproducibility

39. **Lock the conda environment** after install: `mamba env export --no-builds > environment.lock.yml`.
40. **Pin a random seed** in every script that subsamples (HiCRep, APA controls, differential permutations).
41. Use Snakemake `--use-conda` (or `--use-singularity`) so each rule runs in an isolated env.
42. Tag all intermediate files with the assembly version (`hg38_`, `t2t_`, etc.) in the filename.
43. Commit `samples.tsv`, `config/*.yaml`, `environment.yml`, `environment.lock.yml`, and the Snakefile to git. Do not commit results.

## ORACLE integration

44. Always run `09_export_oracle` after loop calling completes. The `.pt` file is the contract with the foundation model.
45. The ORACLE export carries provenance: the `manifest.json` records every input file's SHA-256 hash. Never modify a `.pt` by hand.
46. Sequence-imputed marks must be flagged in the `manifest.json`; the training code can then mask them differently from measured marks.
