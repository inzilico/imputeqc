## Functions to work with vcf files

#' Read vcf files
#'
#' Read vcf and vcf.gz files having both phased or unphased genotypes. The
#' missing genotypes in vcf file can be presented as '.' or './.'
#' \emph{ReadVCF()} is based on \emph{readVcf()} function from
#' VariantAnnotation R package that loads a vcf file as vcf-class object.
#' \emph{ReadVCF()} extracts the genotypes from vcf-class object and converts them
#' in a vector of strings formed by letters A, T, C, G and a symbol of missing data '?'.
#'
#' @param x path/to/filename.vcf, .vcf.gz or .vcf.bgz
#' @param ... an optional parameter \emph{swap} that can be passed to
#'   \emph{Geno2Haps()} function. If \emph{swap} is TRUE, the haplotypes of an
#'   individual are swapped.
#' @return A list with four elements: vcf-class object, a character vector which
#'   elements are strings of letters A, T, C, G, and symbol '?'. The length of
#'   vector is 2N, where N is the number of diploid individuals.
#'   The length of the string equals to the number of markers in vcf file.
#'   Sample ids and marker ids are returned as third and fourth elements
#'   (NULL when not available); they can be passed to
#'   \code{\link{GenerateMaskSet}} to annotate the masks.
#' @export
#'
ReadVCF <- function(x, ...){

  # Check input
  if(!file.exists(x)) stop(x, " doesn't exist!", call. = F)

  # Load vcf
  message("Loading ", x)
  vcf <- VariantAnnotation::readVcf(x)
  # Get genotypes from vcf-class object (matrix variants x samples)
  gt <- VariantAnnotation::geno(vcf)[["GT"]]

  l <- Geno2Haps(gt, ...)

  # Print info about dataset
  message("Haplotypes: ", length(l))
  message("SNPs: ", nchar(l[1]))

  # Make output
  out <- list(vcf = vcf, haps = l,
              samples = colnames(gt), markers = rownames(vcf))
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
#' @param compress logical. If TRUE, the output is written bgzip-compressed
#'   with the .vcf.gz extention. Such files are ready for tabix indexing.
#'   The default is FALSE (plain .vcf as in previous versions).
#' @export
#'
UpdateVCF <- function(haps, pref, n, vcf, compress = FALSE) {

  # Set output filename
  if (compress) {
    fn <- if (is.null(n)) paste0(pref, ".vcf.gz")
          else sprintf("%s.m%s.vcf.gz", pref, n)
  } else {
    fn <- if (is.null(n)) paste0(pref, ".vcf")
          else sprintf("%s.m%s.vcf", pref, n)
  }

  # Remove output file if it exists
  if(file.exists(fn)) file.remove(fn)

  # Convert haplotypes into matrix with genotypes
  gt <- Haps2Geno(haps)

  # Update genotypes of vcf object
  VariantAnnotation::geno(vcf) <- gt

  # Save vcf object; writeVcf does not compress on its own, so bgzip
  # the plain file explicitly when compress = TRUE
  if (compress) {
    tmp <- tempfile(fileext = ".vcf")
    on.exit(unlink(tmp), add = TRUE)
    VariantAnnotation::writeVcf(vcf, tmp)
    Rsamtools::bgzip(tmp, fn, overwrite = TRUE)
  } else {
    VariantAnnotation::writeVcf(vcf, fn)
  }
  message(sprintf("File %s is saved", fn))

}

Haps2Geno <- function(haps){
  # Convert haplotypes to genotypes
  # Args:
  #  haps: vector (or list) with haplotypes
  # Returns:
  #  Matrix with genotypes written as unphased ("/"), variants x samples

  haps <- as.character(unlist(haps, use.names = FALSE))

  # Conver haplotypes into matrix of single alleles (2N x M)
  h <- matrix(unlist(strsplit(haps, "", fixed = TRUE), use.names = FALSE),
              nrow = length(haps), byrow = TRUE)
  h[h == "?"] <- "."

  # Merge haplotypes to get genotypes (N x M)
  ind <- seq(1, nrow(h), 2)
  gt <- matrix(paste0(h[ind, ], "/", h[ind + 1, ]), nrow = nrow(h)/2)

  # Transpose to variants x samples
  t(gt)

}

GetHaps <- function(x) {
  # Read vcf, vcf.gz or vcf.bgz file, convert alleles into letters
  # Args:
  #  x: path/to/filename.{vcf,vcf.gz,vcf.bgz}
  # Returns:
  #  List with haplotypes (alleles given as {A,T,C,G}), sample ids and
  #  marker ids taken from the file

  # Load genotypes as letters (matrix variants x samples)
  gt <- VariantAnnotation::readGT(x, nucleotides = TRUE)

  list(haps = Geno2Haps(gt),
       samples = colnames(gt),
       markers = rownames(gt))
}

Geno2Haps <- function(gt, swap = FALSE) {
  # Convert matrix with genotypes into vector with haplotypes. Symbol of missing
  # values '.' is replaced for '?'
  # Args:
  #  gt: matrix with genotypes (variants x samples)
  #  swap: logic, if TRUE, swap haplotypes of an individual
  # Returns:
  #  Vector with haplotypes

  # Convert genotypes into a vector of haplotypes
  l <- unlist(lapply(seq_len(ncol(gt)), function(j) {
    t <- gsub(".", "?", gt[, j], fixed = TRUE)
    # A genotype written as single '.' becomes a missing diploid genotype
    t[t == "?"] <- "?/?"
    h <- do.call(rbind, strsplit(t, "\\/|\\|"))
    c(paste(h[, 1], collapse = ""), paste(h[, 2], collapse = ""))
  }), use.names = FALSE)

  if(swap) {

    # Swap haplotypes
    message("Swapping haplotypes...")

    m <- matrix(l, ncol = 2, byrow = TRUE)
    l <- as.vector(t(m[, 2:1]))
  }

  return(l)

}
