#!/usr/bin/Rscript
# Script is part of the imputeqc pipeline to estimate number of haplotypes https://github.com/inzilico/imputeqc.
# It reads inp (fastPHASE input file), imputed and masked data and estimates 
# the proportion of wrongly imputed genotypes.

library(imputeqc)

inpuFile <- commandArgs(trailingOnly = T)

# Count errors for one set of test files
errors <- EstimateQuality(origin = inpuFile[1],
                          mask = inpuFile[3], 
                          imputed = inpuFile[2])
