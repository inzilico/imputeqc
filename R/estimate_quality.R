#' Estimate the quality of imputation
#'
#' Counts the discordance of genotypes imputed.
#'
#' A masked genotype is called discordant when its allele pair (order
#' ignored) differs from the allele pair of the original genotype. The
#' discordance is the proportion of such genotypes among all hidden ones.
#' Since version 1.1.0 the comparison is symmetric: a homozygote imputed as
#' heterozygote counts as discordant, as well as the other way round. The
#' allele discordance is also reported: the proportion of mismatched alleles
#' among all alleles of the hidden genotypes. Genotypes that the imputation
#' tool left missing are counted as discordant and are additionally reported
#' in the \emph{n_unfilled} column.
#'
#' The function aligns the imputed and original genotypes by position. For
#' VCF inputs the marker ids (rownames) and sample ids are checked when
#' available: an error is raised when the markers are not identical and in
#' the same order, or when the sample sets differ; a permuted sample order is
#' silently re-ordered. When imputing with a reference panel (e.g. BEAGLE),
#' subset the imputed output back to the markers of the original file first.
#'
#' @param origin path to original unmasked file
#' @param masks path/to/masks.RDS, where \emph{masks.RDS} keeps the masks
#'   created upstream with \emph{\link{GenerateMaskSet}} and saved as
#'   \emph{*.RDS}.
#' @param imputed vector of paths to files imputed. In the case of fastPHASE,
#'   the elements are \emph{path/to/*_genotypes.out}. In the case of vcf files,
#'   the elements are \emph{path/to/*.vcf}, \emph{*.vcf.gz} or \emph{*.vcf.bgz}.
#'   The order of elements of vector and masks should coincide.
#' @param id character or numeric, id of computational experiment. If you run
#'   several calculations with different model parameter to find the optimal
#'   one, you can mark each run with id for the convinience of further
#'   visualization. The argument is optional.
#' @param af_bins optional numeric vector of minor allele frequency breaks.
#'   When given, the discordance is additionally estimated per MAF stratum of
#'   the original data. The breaks must cover the interval from 0 to 0.5, e.g.
#'   \code{c(0, 0.01, 0.05, 0.1, 0.5)}; the strata are intervals open on the
#'   left and closed on the right. The result then contains one row per
#'   imputed file and stratum, the overall row marked as \code{"ALL"} in the
#'   \emph{af_bin} column.
#'
#' @return A data frame with columns \emph{discordance},
#'   \emph{allele_discordance}, \emph{n_masked}, \emph{n_unfilled} and
#'   \emph{mask} (ordinal number of the mask) for every imputed file, and
#'   optional \emph{id}. When \emph{af_bins} is given, the column \emph{af_bin}
#'   is appended and one row per stratum (plus the overall \code{"ALL"} row)
#'   is returned for every imputed file.
#' @export
#'
EstimateQuality <- function(origin, masks, imputed, id = NULL, af_bins = NULL){

  # Load unmasked data set
  g0 <- LoadGenotypes(origin)
  # Dimensions used by the alignment checks
  g0$n_haps <- length(g0$haps)
  g0$m_markers <- nchar(g0$haps[1])

  # Convert to matrix
  c0 <- seq2mat(g0$haps)
  g0$haps <- NULL

  # Output missingness for unmasked data
  ms <- c0 == "?"
  message("Missingness of original data: ", round(sum(ms)/length(ms), 4))
  rm(ms)

  # Load masks
  masks <- readRDS(masks)

  # Minor allele frequency per marker estimated on the original data
  maf <- NULL
  if (!is.null(af_bins)) {
    maf <- MarkerMAF(c0)
  }

  # Initilize output
  out <- vector("list", length(imputed))

  # Loop through all imputed files
  for(i in seq_along(imputed)){

    # Load imputed data and check it is aligned with the original one
    g1 <- LoadGenotypes(imputed[i])
    g1 <- AlignWithOrigin(g0, g1, imputed[i])

    # Convert haplotypes into matrix with alleles
    c1 <- seq2mat(g1$haps)
    g1$haps <- NULL

    # Select mask
    m <- masks[[i]]

    # Convert into TRUE/FALSE and replicate rows
    m <- m == 1
    m <- m[rep(seq_len(nrow(m)), each = 2), , drop = FALSE]

    # Subset original and imputed alleles by mask; each row of the resulting
    # matrices is one hidden genotype (two alleles)
    sel <- which(m)
    a0 <- matrix(c0[sel], ncol = 2, byrow = TRUE)
    a1 <- matrix(c1[sel], ncol = 2, byrow = TRUE)

    # Index of the marker each hidden genotype belongs to (arrayInd columns
    # are row, col; consecutive pairs of selected cells form one genotype)
    geno_marker <- arrayInd(sel, dim(m))[, 2][c(TRUE, FALSE)]

    # Genotype discordance: compare allele pairs with the order ignored
    lo0 <- pmin(a0[, 1], a0[, 2]); hi0 <- pmax(a0[, 1], a0[, 2])
    lo1 <- pmin(a1[, 1], a1[, 2]); hi1 <- pmax(a1[, 1], a1[, 2])
    gt_disc <- (lo0 != lo1) | (hi0 != hi1)

    # Allele discordance
    al_disc <- (lo0 != lo1) + (hi0 != hi1)

    # Genotypes the imputation tool left missing
    unfilled <- (a1[, 1] == "?") | (a1[, 2] == "?")

    n <- nrow(a0)
    disc <- round(sum(gt_disc)/n, 6)
    message("Discordance: ", disc)

    # Make the output
    row_i <- data.frame(discordance = disc,
                        allele_discordance = round(sum(al_disc)/(2*n), 6),
                        n_masked = n,
                        n_unfilled = sum(unfilled),
                        mask = i,
                        stringsAsFactors = FALSE)

    if (!is.null(af_bins)) {
      bin <- cut(maf[geno_marker], breaks = af_bins, right = TRUE,
                 include.lowest = TRUE)
      if (any(is.na(bin)))
        warning(sprintf(paste("%s hidden genotypes fall outside af_bins and",
                              "are excluded from the per-stratum rows"),
                        sum(is.na(bin))),
                call. = FALSE, immediate. = TRUE)
      by_bin <- do.call(rbind, lapply(levels(bin), function(b) {
        k <- which(bin == b)
        data.frame(discordance = if (length(k)) round(sum(gt_disc[k])/length(k), 6) else NA_real_,
                   allele_discordance = if (length(k)) round(sum(al_disc[k])/(2*length(k)), 6) else NA_real_,
                   n_masked = length(k),
                   n_unfilled = if (length(k)) sum(unfilled[k]) else 0L,
                   mask = i,
                   af_bin = b,
                   stringsAsFactors = FALSE)
      }))
      row_i <- rbind(cbind(row_i, af_bin = "ALL"), by_bin)
    }

    if (!is.null(id)) row_i$id <- id
    out[[i]] <- row_i
  }

  do.call(rbind, out)
}

MarkerMAF <- function(alleles) {
  # Minor allele frequency per marker from a matrix of single-character
  # alleles (2N x M). Works with any allele coding (letters or numbers) as
  # the files keep: the frequency of the non-major characters is returned,
  # which equals the minor allele frequency for biallelic markers.
  vapply(seq_len(ncol(alleles)), function(j) {
    v <- alleles[, j]
    v <- v[v != "?"]
    if (!length(v)) return(NA_real_)
    1 - max(table(v))/length(v)
  }, numeric(1))
}

AlignWithOrigin <- function(g0, g1, label) {
  # Checks that imputed genotypes can be compared with the original ones
  # Args:
  #   g0, g1: lists returned by LoadGenotypes for the original and imputed
  #     files; g1 is returned re-ordered to the sample order of g0 when needed
  #   label: path of the imputed file, for error messages

  if (length(g1$haps) != g0$n_haps ||
      any(nchar(g1$haps) != g0$m_markers))
    stop(sprintf(paste("Imputed file %s has %s haplotypes of %s markers,",
                       "while the original file has %s of %s. Align the files",
                       "first: with reference-based imputation subset the",
                       "imputed output back to the markers of the original",
                       "file."),
                 label, length(g1$haps), nchar(g1$haps[1]),
                 g0$n_haps, g0$m_markers),
         call. = FALSE)

  if (!is.null(g0$samples) && !is.null(g1$samples)) {
    if (setequal(g0$samples, g1$samples)) {
      if (!identical(g0$samples, g1$samples)) {
        message("Re-ordering samples of ", label, " to match the original file")
        ord <- match(g0$samples, g1$samples)
        m <- matrix(g1$haps, ncol = 2, byrow = TRUE)
        g1$haps <- as.vector(t(m[ord, , drop = FALSE]))
        g1$samples <- g1$samples[ord]
      }
    } else {
      stop(sprintf("Sample sets differ between the original file and %s: %s",
                   label,
                   paste(setdiff(union(g0$samples, g1$samples),
                                 intersect(g0$samples, g1$samples)),
                         collapse = ", ")),
           call. = FALSE)
    }
  }

  if (!is.null(g0$markers) && !is.null(g1$markers)) {
    bad <- which(g0$markers != g1$markers)
    if (length(bad) > 0)
      stop(sprintf(paste("Markers of %s do not coincide with the original",
                         "file: first mismatch at position %s (%s vs %s).",
                         "The markers must be identical and in the same order;",
                         "subset the imputed file accordingly."),
                   label, bad[1],
                   g0$markers[bad[1]], g1$markers[bad[1]]),
           call. = FALSE)
  }

  g1
}

# suppress R CMD check note about the ggplot non-standard evaluation
utils::globalVariables("discordance")

#' Boxpolt discordance
#'
#' Draw boxplots of discordance estimated for different model parameters.
#'
#' @param fname path to dataframe with columns 'discordance' and 'id'
#' @param tl optional title of the plot
#' @param stl optional subtitle of the plot
#' @param id optional x-axis label
#'
#' @return ggplot object
#' @import ggplot2
#' @export
#'
PlotDiscordance <- function(fname, tl = NULL, stl = NULL, id = NULL){

  if(!file.exists(fname)) stop((sprintf("File %s doesn't exist", fname)))

  # Load data set
  df <- utils::read.table(fname, header = T)

  # Keep the overall rows when the table comes from af-stratified estimation
  if ("af_bin" %in% colnames(df)) {
  df <- df[df$af_bin == "ALL" | is.na(df$af_bin), , drop = FALSE]
  }

  df$id <- factor(df$id)

  # Plot
  p <- ggplot(df, aes(y = discordance, x = id)) +
    geom_boxplot(varwidth = TRUE) +
    labs(title = ifelse(is.null(tl),"",tl), subtitle = ifelse(is.null(stl), "", stl),
         y = "Discordance", x = ifelse(is.null(id), "", id)) +
    scale_x_discrete(limits = levels(df$id)) +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
          title = element_text(size = 16))

  # Return ggplot object
  return(p)

}

seq2mat <- function(vec){
  # Converts sequences into matrix (one row per sequence, one column per allele)

  matrix(unlist(strsplit(vec, "", fixed = TRUE), use.names = FALSE),
         nrow = length(vec), byrow = TRUE)

}

LoadGenotypes <- function(input){
  # Load data and convert them into haplotypes. fastPHASE as well as vcf
  # (plain, gzipped or bgzipped) formats are accessable. Sample and marker ids
  # are extracted from VCF files when present.
  # Args:
  #  input: path/to/filename.{inp,out,vcf,vcf.gz,vcf.bgz}
  # Returns:
  #  List with haplotypes, sample ids and marker ids (the last two are NULL
  #  for fastPHASE files)

  # Get extention of input file
  m <- regexpr("\\.([[:alnum:]]+)$", input)
  ext <- tolower(regmatches(input, m))

  g <- switch(ext,
              .vcf = GetHaps(input),
              .gz = GetHaps(input),
              .bgz = GetHaps(input),
              .inp = list(haps = ReadFastPHASE(input), samples = NULL, markers = NULL),
              .out = list(haps = ReadFastPHASE(input), samples = NULL, markers = NULL))

  if(is.null(g)) stop("Haplotypes aren't loaded! Unsupported file type: ", input,
                      call. = F)

  return(g)

}
