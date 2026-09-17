# Plain-R test suite for imputeqc (runs under R CMD check without testthat).
# Every assertion is a stopifnot(); the script prints "imputeqc tests: OK".

suppressMessages(library(imputeqc))
set.seed(42)

td <- file.path(tempdir(), "imputeqc_tests")
dir.create(td, recursive = TRUE, showWarnings = FALSE)

pass <- function(label) cat("PASS:", label, "\n")

# --- fixture: a small VCF with known genotypes -------------------------------
# 6 samples x 10 markers, alt counts 1..12 over 12 alleles -> varied MAF
samples <- sprintf("S%d", 1:6)
markers <- sprintf("rs%d", 1:10)
gts <- matrix(c(
  "0/0","0/0","0/0","0/0","0/0","0/1",  # rs1
  "0/0","0/0","0/1","0/0","0/1","0/0",  # rs2
  "0/0","0/0","0/1","0/1","0/0","0/0",  # rs3
  "0/0","0/1","0/1","0/0","0/1","0/1",  # rs4
  "0/0","0/1","0/1","0/1","0/1","0/0",  # rs5
  "0/1","0/1","0/1","0/1","0/1","0/1",  # rs6
  "1/1","0/0","0/0","0/0","0/0","0/0",  # rs7
  "1/1","0/1","0/0","0/0","0/0","0/0",  # rs8
  "1/1","1/1","0/1","0/0","0/0","0/0",  # rs9
  "1/1","1/1","1/1","0/1","0/0","0/0"), # rs10
  ncol = 10, byrow = FALSE)
stopifnot(all(dim(gts) == c(6, 10)))
gts[3, 8] <- "./."   # S3 rs8 missing
gts[5, 9] <- "./."   # S5 rs9 missing

write_test_vcf <- function(file, g, samples, markers,
                           shuffle_samples = FALSE, rename_marker = NULL,
                           drop_sample = FALSE, rename_sample = FALSE) {
  j <- seq_along(samples)
  if (shuffle_samples) j <- rev(j)
  if (drop_sample) j <- j[-1]
  out_samples <- samples[j]
  if (rename_sample) out_samples[length(out_samples)] <- "S9"
  mk <- markers
  if (!is.null(rename_marker)) mk[rename_marker] <- "rsWRONG"
  hdr <- c(
    "##fileformat=VCFv4.2",
    "##contig=<ID=1,length=100000>",
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
  body <- vapply(seq_along(mk), function(i) {
    paste(c(paste("1", 100 + i, mk[i], if (i %% 2) "A" else "C",
                 if (i %% 2) "G" else "T", ".", "PASS", ".", "GT", sep = "\t"),
            g[j, i]), collapse = "\t")
  }, character(1))
  writeLines(c(hdr, paste(c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER",
                            "INFO","FORMAT", out_samples), collapse = "\t"),
               body), file)
  invisible(file)
}

origin_vcf <- file.path(td, "origin.vcf")
write_test_vcf(origin_vcf, gts, samples, markers)

# --- ReadVCF ------------------------------------------------------------------
x <- ReadVCF(origin_vcf)
stopifnot(length(x$haps) == 12)                # 2N haplotypes
stopifnot(nchar(x$haps[1]) == 10)              # of M markers
stopifnot(identical(x$samples, samples))
stopifnot(identical(x$markers, markers))
# the two missing genotypes show as '?' at both haplotypes of the carrier;
# this code path keeps the allele symbols of the GT field (0/1 here)
stopifnot(substring(x$haps[2 * which(samples == "S3") - 1], 8, 8) == "?",
          substring(x$haps[2 * which(samples == "S3")], 8, 8) == "?",
          substring(x$haps[2 * which(samples == "S5") - 1], 9, 9) == "?",
          substring(x$haps[2 * which(samples == "S5")], 9, 9) == "?")
stopifnot(sum(vapply(x$haps, function(s) sum(strsplit(s, "")[[1]] == "?"), 0)) == 4)
# swap = TRUE exchanges the two haplotypes of every individual
xs <- ReadVCF(origin_vcf, swap = TRUE)
h <- matrix(x$haps, ncol = 2, byrow = TRUE)
hs <- matrix(xs$haps, ncol = 2, byrow = TRUE)
stopifnot(identical(h[, 1], hs[, 2]), identical(h[, 2], hs[, 1]))
pass("ReadVCF")

# --- GenerateMaskSet ----------------------------------------------------------
masks <- GenerateMaskSet(x$haps, n = 3, p = 0.4, samples = samples, markers = markers)
stopifnot(length(masks) == 3)
stopifnot(all(vapply(masks, function(m) all(dim(m) == c(6, 10)), logical(1))))
stopifnot(identical(rownames(masks[[1]]), samples),
          identical(colnames(masks[[1]]), markers))
# masks within the set are disjoint and avoid originally missing genotypes
tot <- masks[[1]] + masks[[2]] + masks[[3]]
stopifnot(all(tot <= 1), tot[3, 8] == 0, tot[5, 9] == 0)
# each mask hides about p of the available genotypes
n_avail <- 6 * 10 - 2
stopifnot(all(vapply(masks, function(m) abs(sum(m)/n_avail - 0.4) < 0.1,
                     logical(1))))
# wrong annotation lengths are rejected
r <- try(GenerateMaskSet(x$haps, 2, 0.1, samples = samples[1:4]), silent = TRUE)
stopifnot(inherits(r, "try-error"))
# marker hiding masks whole markers
mm <- GenerateMaskSet(x$haps, n = 2, p = 0.2, type = "marker")
stopifnot(all(mm[[1]] %in% c(0, 1)),
          all(colSums(mm[[1]]) %in% c(0, 6)),
          sum(colSums(mm[[1]]) > 0) == 2)
pass("GenerateMaskSet")

# --- ApplyMasks + UpdateVCF ---------------------------------------------------
saveRDS(masks, file.path(td, "masks.RDS"))
ApplyMasks(x$haps, masks, file.path(td, "masked"), vcf = x$vcf)
stopifnot(all(file.exists(file.path(td, sprintf("masked.m%d.vcf", 1:3)))))
xm <- ReadVCF(file.path(td, "masked.m1.vcf"))
stopifnot(identical(xm$samples, samples), identical(xm$markers, markers))
# hidden genotypes became missing; all other alleles are untouched
allele_mat <- function(haps) matrix(unlist(strsplit(haps, "", fixed = TRUE),
                                           use.names = FALSE),
                                    nrow = length(haps), byrow = TRUE)
A0 <- allele_mat(x$haps)
A1 <- allele_mat(xm$haps)
m1r <- (masks[[1]] == 1)[rep(seq_len(6), each = 2), , drop = FALSE]
stopifnot(all(A1[m1r] == "?"))
stopifnot(all(A0[!m1r] == A1[!m1r]))
# bgzip-compressed output is a real gzip file with the same content
UpdateVCF(x$haps, file.path(td, "comp"), 1, x$vcf, compress = TRUE)
magic <- readBin(file.path(td, "comp.m1.vcf.gz"), "raw", n = 2)
stopifnot(identical(as.character(magic), c("1f", "8b")))
xc <- ReadVCF(file.path(td, "comp.m1.vcf.gz"))
stopifnot(identical(xc$haps, x$haps))
pass("ApplyMasks + UpdateVCF")

# --- EstimateQuality on a simulated imputation -------------------------------
# the imputed file carries the true genotypes except: one hom imputed as het,
# one het imputed as hom, and one genotype the imputer left unfilled
cells <- which(masks[[1]] == 1, arr.ind = TRUE)
# pairwise (row, col) indexing: gts[rvec, cvec] would give the cartesian product
truths <- gts[cbind(cells[, "row"], cells[, "col"])]
hom_row <- which(truths == "0/0")[1]
het_row <- which(truths == "0/1")[1]
# any remaining masked cell can play the unfilled one (the masks never hide
# the originally missing genotypes)
others <- setdiff(seq_len(nrow(cells)), c(hom_row, het_row))
unf_row <- others[1]
stopifnot(truths[unf_row] != "./.")
hom_cell <- cells[hom_row, , drop = TRUE]
het_cell <- cells[het_row, , drop = TRUE]
unf_cell <- cells[unf_row, , drop = TRUE]
g_imp <- gts
g_imp[hom_cell[1], hom_cell[2]] <- "0/1"
g_imp[het_cell[1], het_cell[2]] <- "1/1"
g_imp[unf_cell[1], unf_cell[2]] <- "./."
imp_vcf <- file.path(td, "imputed.m1.vcf")
write_test_vcf(imp_vcf, g_imp, samples, markers)

eq <- EstimateQuality(origin = origin_vcf, masks = file.path(td, "masks.RDS"),
                      imputed = imp_vcf)
stopifnot(nrow(eq) == 1, eq$mask == 1)
stopifnot(eq$n_masked == sum(masks[[1]]), eq$n_unfilled == 1)
stopifnot(eq$discordance == round(3/eq$n_masked, 6))   # 2 wrong + 1 unfilled
# allele errors: hom->het 1, het->hom 1, unfilled 2
stopifnot(eq$allele_discordance == round(4/(2 * eq$n_masked), 6))
# the metric is symmetric: a homozygote imputed as het counts as discordant
stopifnot(gts[hom_cell[1], hom_cell[2]] == "0/0")
pass("EstimateQuality (simulated imputation)")

# --- allele-frequency strata --------------------------------------------------
eqb <- EstimateQuality(origin = origin_vcf, masks = file.path(td, "masks.RDS"),
                       imputed = imp_vcf, af_bins = c(0, 0.25, 0.5))
stopifnot(nrow(eqb) == 3)                                # ALL + two strata
all_row <- eqb[eqb$af_bin == "ALL", , drop = FALSE]
stopifnot(all_row$discordance == eq$discordance,
          all_row$n_masked == eq$n_masked,
          sum(eqb$n_masked[eqb$af_bin != "ALL"]) == eq$n_masked)
pass("EstimateQuality (af strata)")

# --- validation of alignment --------------------------------------------------
# a file with an extra marker must be rejected
g_extra <- cbind(gts, gts[, 1])
extra_vcf <- file.path(td, "imputed_extra.vcf")
mk_extra <- c(markers, "rs11")
write_test_vcf(extra_vcf, g_extra, samples, mk_extra)
r <- try(EstimateQuality(origin_vcf, file.path(td, "masks.RDS"), extra_vcf),
         silent = TRUE)
stopifnot(inherits(r, "try-error"),
          grepl("Align the files", r, fixed = TRUE))
# a renamed marker must be rejected with the mismatch reported
ren_vcf <- file.path(td, "imputed_renamed.vcf")
write_test_vcf(ren_vcf, gts, samples, markers, rename_marker = 5)
r <- try(EstimateQuality(origin_vcf, file.path(td, "masks.RDS"), ren_vcf),
         silent = TRUE)
stopifnot(inherits(r, "try-error"), grepl("rsWRONG", r, fixed = TRUE))
# a dropped sample changes the dimensions and is rejected by the size check
drop_vcf <- file.path(td, "imputed_dropped.vcf")
write_test_vcf(drop_vcf, gts, samples, markers, drop_sample = TRUE)
r <- try(EstimateQuality(origin_vcf, file.path(td, "masks.RDS"), drop_vcf),
         silent = TRUE)
stopifnot(inherits(r, "try-error"), grepl("Align the files", r, fixed = TRUE))
# a foreign sample id at the same count is rejected by the sample check
foreign_vcf <- file.path(td, "imputed_foreign.vcf")
write_test_vcf(foreign_vcf, gts, samples, markers, rename_sample = TRUE)
r <- try(EstimateQuality(origin_vcf, file.path(td, "masks.RDS"), foreign_vcf),
         silent = TRUE)
stopifnot(inherits(r, "try-error"), grepl("Sample sets differ", r, fixed = TRUE))
# a permuted sample order is re-aligned by id: corrupt one masked genotype
# and check the discordance is exactly 1/n
mc <- cells[1, , drop = TRUE]
g_perm <- gts
g_perm[mc[1], mc[2]] <- if (gts[mc[1], mc[2]] == "0/1") "1/1" else "0/1"
perm_vcf <- file.path(td, "imputed_permuted.vcf")
write_test_vcf(perm_vcf, g_perm, samples, markers, shuffle_samples = TRUE)
eqp <- EstimateQuality(origin_vcf, file.path(td, "masks.RDS"), perm_vcf)
stopifnot(eqp$discordance == round(1/eqp$n_masked, 6), eqp$n_unfilled == 0)
pass("alignment checks")

# --- fastPHASE roundtrip ------------------------------------------------------
inp <- file.path(td, "origin")            # WriteFastPHASE adds the .inp suffix
WriteFastPHASE(x$haps, inp)
inp <- paste0(inp, ".inp")
stopifnot(identical(ReadFastPHASE(inp), x$haps))
eqf <- EstimateQuality(origin = inp, masks = file.path(td, "masks.RDS"),
                       imputed = inp)
stopifnot(eqf$discordance == 0, eqf$n_unfilled == 0)
pass("fastPHASE roundtrip")

# --- WriteMaskSet -------------------------------------------------------------
mask_tsv <- file.path(td, "masks.tsv")
WriteMaskSet(masks, mask_tsv)
tab <- read.table(mask_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
stopifnot(identical(colnames(tab), c("mask", "sample", "marker")),
          nrow(tab) == sum(vapply(masks, sum, numeric(1))))
k2 <- tab[tab$mask == 2, ]
cells2 <- which(masks[[2]] == 1, arr.ind = TRUE)
stopifnot(identical(sort(paste(k2$sample, k2$marker, sep = "|")),
                    sort(paste(samples[cells2[, "row"]],
                               markers[cells2[, "col"]], sep = "|"))))
pass("WriteMaskSet")

# --- PlotDiscordance ----------------------------------------------------------
eq5 <- EstimateQuality(origin_vcf, file.path(td, "masks.RDS"),
                       rep(imp_vcf, 2), id = 5)
eq10 <- EstimateQuality(origin_vcf, file.path(td, "masks.RDS"),
                        rep(imp_vcf, 2), id = 10)
tab2 <- rbind(eq5, eq10)
stopifnot(nrow(tab2) == 4)                      # 2 ids x 2 masks
plot_file <- file.path(td, "discordance.tsv")
write.table(tab2, plot_file, sep = "\t", row.names = FALSE)
p <- PlotDiscordance(plot_file)
stopifnot(inherits(p, "ggplot"), identical(levels(p$data$id), c("5", "10")))
# af-stratified tables are reduced to their ALL rows by the plot
eqb_id <- EstimateQuality(origin_vcf, file.path(td, "masks.RDS"), imp_vcf,
                          id = 5, af_bins = c(0, 0.25, 0.5))
plot_file_b <- file.path(td, "discordance_bins.tsv")
write.table(eqb_id, plot_file_b, sep = "\t", row.names = FALSE)
pb <- PlotDiscordance(plot_file_b)
stopifnot(nrow(pb$data) == 1)
pass("PlotDiscordance")

cat("imputeqc tests: OK\n")
