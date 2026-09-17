What is it
----------

*imputeqc* is an R package and accompanied scripts to estimate the quality of imputation of genotypes that was made with [fastPHASE](http://scheet.org/software.html) and [BEAGLE](https://faculty.washington.edu/browning/beagle/beagle.html) softwares. Any other tools that support *.inp fastPHASE or VCF files can be also considered. The package is based on masked data analysis. 

Possible applications
---------------------

1. Estimation of the error of gynotype imputation.

2. Optimization of the imputation model parameters, e.g., the number of haplotype clusters. The parameter can be further used for the search of signatures of selection with [hapFLK](https://forge-dga.jouy.inra.fr/projects/hapflk) test.

3. Testing different reference panels for imputation.

4. Benchmarking of different imputation softwares and strategies.

Ho to install from GitHub
--------------------------

Run from R.

1. Make sure you have [devtools](https://github.com/r-lib/devtools) package installed. 

```
install.packages("devtools")
```

2. Install dependencies

```
install.packages("BiocManager")
BiocManager::install("VariantAnnotation")
```

3. Install *imputeqc*

```
devtools::install_github("inzilico/imputeqc", build_vignettes = TRUE)
```

How to use
----------

Read a vignette [How to Select the Number of Clusters for fastPHASE](https://htmlpreview.github.io/?https://github.com/inzilico/imputeqc/blob/master/vignettes/k_selection.html). 

* On a local machine, the vignette can be accessed as follow: 
```
browseVignettes("imputeqc")
```    
* On remote machine, the vignette can be opened in the "Help" tab of RStudio:
```
vignette("k_selection")
```

Quality of BEAGLE imputation (VCF files)
----------------------------------------

The package works with VCF files, which is convenient for BEAGLE. The
typical workflow:

```r
library(imputeqc)

# 1. Load the genotypes and hide 1% of them in each of the 3 masks
x <- ReadVCF("genotypes.vcf.gz")
set.seed(42)
masks <- GenerateMaskSet(x$haps, n = 3, p = 0.01,
                         samples = x$samples, markers = x$markers)
saveRDS(masks, "masks.RDS")
ApplyMasks(x$haps, masks, pref = "masked", vcf = x$vcf)
# optional: export the hidden genotypes for use outside R
WriteMaskSet(masks, "masks.tsv")

# 2. Impute every masked file with BEAGLE (run in shell):
#    java -jar beagle.jar gt=masked.m1.vcf ref=reference.vcf.gz out=imputed.m1
#    ...
#    With reference-based imputation subset every output back to the
#    markers of the original file (e.g. bcftools view -T) before step 3.

# 3. Estimate the discordance of the imputed genotypes
eq <- EstimateQuality(origin = "genotypes.vcf.gz",
                      masks = "masks.RDS",
                      imputed = c("imputed.m1.vcf.gz",
                                  "imputed.m2.vcf.gz",
                                  "imputed.m3.vcf.gz"))
# optionally per minor allele frequency stratum
eqb <- EstimateQuality(origin = "genotypes.vcf.gz",
                       masks = "masks.RDS",
                       imputed = c("imputed.m1.vcf.gz",
                                   "imputed.m2.vcf.gz",
                                   "imputed.m3.vcf.gz"),
                       af_bins = c(0, 0.01, 0.05, 0.1, 0.5))
```

See the vignette `beagle_quality` for a complete runnable example.

Changes in version 1.1.0 are listed in [NEWS.md](NEWS.md).

License
-------
[MIT](https://en.wikipedia.org/wiki/MIT_License)

Citing
------
Khvorykh GV, Khrunin AV. imputeqc: an R package for assessing imputation quality of genotypes and optimizing imputation parameters. BMC Bioinformatics. 2020;21(Suppl 12):304. Published 2020 Jul 24. doi:10.1186/s12859-020-03589-0. [\[pubmed\]](https://pubmed.ncbi.nlm.nih.gov/32703240/), [\[pdf\]](https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/s12859-020-03589-0)

Contacts
--------
Gennady Khvorykh, a bioinformatician, [inzilico.com](http://inzilico.com)

Interested in contributing to the project? Suggestions, questions, and comments are open! Feel free [to drop me the message](http://www.inzilico.com/contacts/).
