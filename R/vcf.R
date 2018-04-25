## Functions to work with vcf files

#' Read vcf files
#'
#' Read vcf and vcf.gz files having both phased or unphased genotypes. The
#' missing genotypes can be presented as '.' or './.'
#'
#' @param x path/to/filename.vcf or vcf.gz
#' @param ... an optional parameter \emph{swap} that can be passed to
#'   \emph{Geno2Haps()} function. If \emph{swap} is TRUE, the haplotypes of an
#'   individual are swapped.
#' @return A list with two elements: vcf-class object and a character vector which
#'   elements are strings of letters representing haplotypes. The length of
#'   vector is 2N, where N is the number of diploid individuals.
#' @export
#'
ReadVCF <- function(x, ...){

  # Check input
  if(!file.exists(x)) stop(x, " doesn't exist!", call. = F)

  # Load vcf
  message("Loading ", x)
  vcf <- VariantAnnotation::readVcf(x)
  # Get genotypes from vcf-class object
  gt <- VariantAnnotation::geno(vcf)@listData$GT

  l <- Geno2Haps(gt, ...)

  # Print info about dataset
  message("Haplotypes: ", length(l))
  message("SNPs: ", nchar(l[1]))

  # Make output
  out <- list(vcf, l)
  names(out) <- c("vcf", "haps")

  return(out)

}

#' Update genotypes of vcf file
#'
#' Update GT fields of vcf file by haplotypes provided. Updated file is saved
#' under the name, which can be numbered to distinguish files in a set. It is
#' useful, when you create a set of vcf files where GT fields are masked in
#' different ways. These files will be further imputed. The updated genotypes
#' are written as unphased.
#'
#' @param haps character vector with haplotypes
#' @param pref path/to/pref, where `pref` is a basename of output file. vcf
#'   extention is added.
#' @param n ordinal number to mark updated files automatically. If provided,
#'   \emph{m} followed by this number is added to the basename of the file.
#' @param vcf vcf-class object loaded with \emph{readVcf()} function from
#'   VariantAnnotation package. The genotypes of this object we want to update.
#' @export
#'
UpdateVCF <- function(haps, pref, n, vcf) {

  # Set output filename
  fn <- ifelse(is.null(n), paste0(pref, ".vcf"), sprintf("%s.m%s.vcf", pref, n))

  # Remove output file if it exists
  if(file.exists(fn)) file.remove(fn)

  # Convert haplotypes into matrix with genotypes
  gt <- Haps2Geno(haps)

  # Update genotypes of vcf object
  VariantAnnotation::geno(vcf) <- gt

  # Save vcf object
  VariantAnnotation::writeVcf(vcf, fn)
  message(sprintf("File %s is saved", fn))

}

Haps2Geno <- function(haps){
  # Convert haplotypes to genotypes
  # Args:
  #  haps: vector with haplotypes
  # Returns:
  #  Matrix with genotypes written as unphased ("/")

  # Conver haplotypes into genotypes
  h <- plyr::laply(haps, function(x) {
    t <- unlist(strsplit(x, split = ""))
    t <- gsub("\\?", ".", t)
    t
    })

  # Get odd rows
  ind <- seq(1, nrow(h), 2)

  # Merge haplotypes to get genotypes
  gt <- paste(h[ind, ], h[-ind, ], sep = "/")
  gt <- matrix(gt, nrow = nrow(h)/2)

  # Transpose
  gt <- t(gt)

  return(gt)
}

GetHaps <- function(x) {
  # Read vcf and vcf.gz files, convert alleles into letters
  # Args:
  #  x: path/to/filename.{vcf,vcf.gz)
  # Returns:
  #  Vector with haplotypes, where alleles are given as {A,T,C,G}

  # Load genotypes as letters
  gt <- VariantAnnotation::readGT(x, nucleotides = T)

  # Convert genotypes into haplotypes
  Geno2Haps(gt)

}

Geno2Haps <- function(gt, swap = FALSE) {
  # Convert matrix with genotypes into vector with haplotypes. Symbol of missing
  # values '.' is replaced for '?'
  # Args:
  #  gt: matrix with genotypes
  #  swap: logic, if TRUE, swap haplotypes of an individual
  # Returns:
  #  Vector with haplotypes

  # Convert genotypes into a vector of haplotypes.
  l <- plyr::alply(gt, 2, function(x) {
    t <- gsub(pattern = "^\\.$", replacement = "?/?", x)
    t <- gsub(pattern = "\\.", replacement = "?", x)
    t <- strsplit(t, "\\/|\\|")
    h <- do.call(rbind, t)
    plyr::alply(h, 2, paste0, collapse = "")

  })

  l <- unlist(l)
  names(l) <- NULL

  if(swap) {

    # Swap haplotypes
    message("Swapping haplotypes...")

    nc <- length(l)

    s <- vector()
    s1 <- seq(1, nc, 2)
    s2 <- seq(2, nc, 2)

    for (i in seq_along(s1)) {
      s <- append(s, s2[i])
      s <- append(s, s1[i])
    }

    l <- l[s]
  }

  return(l)

}

