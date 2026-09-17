# imputeqc 1.1.0

## Breaking changes

- `EstimateQuality` uses a **symmetric genotype comparison**: a genotype is
  discordant when its allele pair (order ignored) differs from the original
  one. Previously a homozygous genotype imputed as heterozygous was counted
  as concordant, which understated the error by a homozygosity-dependent
  amount. Discordance values from this version are therefore higher and not
  directly comparable with <= 1.0 releases.
- The data frame returned by `EstimateQuality` now has the columns
  `discordance`, `allele_discordance`, `n_masked`, `n_unfilled`, `mask`
  (and optional `id`, `af_bin`) instead of the former `discordance` +
  optional `id`.

## New features

- `EstimateQuality` gained the `af_bins` argument for minor-allele-frequency
  stratified estimation of the discordance.
- `GenerateMaskSet` accepts `samples` and `markers` to annotate the masks
  with sample and marker ids (dimnames).
- New function `WriteMaskSet` exports the hidden genotypes of all masks as
  a tab-separated table (`mask`, `sample`, `marker`), making the masks
  portable outside R.
- `EstimateQuality` validates the alignment of the imputed and original
  files: dimensions are always checked; for VCF inputs sample ids and marker
  ids are verified when available (a permuted sample order is re-aligned,
  mismatched markers or sample sets stop with an informative error).
- `UpdateVCF` gained the `compress` argument to write bgzip-compressed
  `.vcf.gz` output (ready for tabix).
- `.vcf.bgz` files are accepted on input.

## Bug fixes

- `Geno2Haps`: a genotype written as a single `.` is now treated as a fully
  missing diploid genotype (`?/?`); previously the first `gsub` result was
  overwritten and such genotypes broke the haplotype assembly.
- `PlotDiscordance`: the x axis levels are inferred from the data instead
  of the hardcoded fastPHASE K values; af-stratified tables are reduced to
  their overall (`ALL`) rows.

## Internal

- The `plyr` dependency is dropped; all operations are base R and the
  hot loops (mask application, haplotype conversion, genotype comparison)
  are vectorized, which speeds up large files considerably.
- A plain-R test suite (`tests/imputeqc_tests.R`) covers reading, masking,
  the discordance metrics, alignment checks, the fastPHASE roundtrip, mask
  export and plotting.
- New vignette `beagle_quality`: masked-data analysis with BEAGLE VCF files.
