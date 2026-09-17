#!/usr/bin/env Rscript
# Generates the example files of the beagle_quality vignette:
#   beagle.vcf                  20 samples x 200 SNPs, varied MAF
#   beagle_imputed.m{1,2,3}.vcf "imputations" of the masks: true genotypes
#                               with 2% of the hidden genotypes corrupted
#                               and one genotype left unfilled
suppressMessages(library(imputeqc))

set.seed(20260917)
n_samples <- 20
n_markers <- 200

samples <- sprintf("ID%02d", seq_len(n_samples))
markers <- sprintf("rsV%03d", seq_len(n_markers))
refs <- sample(c("A", "C"), n_markers, replace = TRUE)
alts <- ifelse(refs == "A", "G", "T")
# alt allele frequencies spread from ~0.005 to 0.5
p <- seq(0.005, 0.5, length.out = n_markers)

gt <- matrix("0/0", nrow = n_samples, ncol = n_markers)
for (j in seq_len(n_markers)) {
  a <- rbinom(n_samples, 2, p[j])
  gt[, j] <- c("0/0", "0/1", "1/1")[a + 1]
}
gt[sample(n_samples * n_markers, 8)] <- "./."   # a little missingness

vcf_lines <- c(
  "##fileformat=VCFv4.2",
  "##contig=<ID=1,length=1000000>",
  '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
  paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
          "FORMAT", samples), collapse = "\t"),
  vapply(seq_len(n_markers), function(j)
    paste(c(paste("1", 1000 + j, markers[j], refs[j], alts[j], ".", "PASS",
                  ".", "GT", sep = "\t"), gt[, j]), collapse = "\t"),
    character(1)))

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
writeLines(vcf_lines, "inst/extdata/beagle.vcf")

# --- the "imputations" --------------------------------------------------------
x <- ReadVCF("inst/extdata/beagle.vcf")
set.seed(42)
masks <- GenerateMaskSet(x$haps, n = 3, p = 0.05,
                         samples = samples, markers = markers)

for (k in seq_along(masks)) {
  h <- matrix(x$haps, ncol = 2, byrow = TRUE)   # truth, 0/1 coded
  cells <- which(masks[[k]] == 1, arr.ind = TRUE)
  # corrupt 2% of the hidden genotypes by flipping the first allele
  flip <- cells[sample(nrow(cells), round(0.02 * nrow(cells))), , drop = FALSE]
  for (r in seq_len(nrow(flip))) {
    i <- flip[r, "row"]; j <- flip[r, "col"]
    ch <- strsplit(h[i, 1], "", fixed = TRUE)[[1]]
    ch[j] <- if (ch[j] == "0") "1" else "0"
    h[i, 1] <- paste(ch, collapse = "")
  }
  # leave one hidden genotype unfilled
  i0 <- cells[1, "row"]; j0 <- cells[1, "col"]
  for (s in 1:2) {
    ch <- strsplit(h[i0, s], "", fixed = TRUE)[[1]]
    ch[j0] <- "?"
    h[i0, s] <- paste(ch, collapse = "")
  }
  haps <- as.vector(t(h))
  UpdateVCF(haps, "inst/extdata/beagle_imputed", k, x$vcf)
}
