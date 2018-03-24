#' Estimate the quality of imputation
#'
#' \code{EstimateQuality} estimates the error of allele and genotype imputation.
#'
#' @param origin path/to/filename.inp, where \code{filename.inp} is the original
#'   (unmasked) fastPHASE file
#' @param masks path/to/masks.RDS, where \code{masks.RDS} keeps the masks
#'   created upstream with \code{\link{GenerateMaskSet}} and saved as
#'   \code{*.RDS}.
#' @param imputed vector of \code{path/to/*_genotypes.out}. The files have
#'   genotypes imputed with fastPHASE. They should be in the same order as masks
#'   were generated.
#' @param id character, id of computational experiment. In case you run several
#'   calculations with different model parameter to find the best one, you can
#'   mark each run with id for the convinience of further visualization. The
#'   argument is optional.
#'
#' @return A data frame with three columns: "alleles", "genotypes", and "id" if
#'   provided. The first two contains the errors, the third one - id of the
#'   experiment. The error is counted as \code{(1 - accuracy)}, where
#'   \code{accuracy} is a proportion of correctly imputed genotypes and alleles.
#'   By correctly imputed genotype we mean that both alleles coincided with the
#'   original ones. The function returns the values for one set of test files.
#' @export
#'
#' @examples
EstimateQuality <- function(origin, masks, imputed, id = NULL){

  # Load original data set
  g0 <- ReadFastPHASE(origin)

  # Convert to matrix
  c0 <- seq2mat(g0)
  rm(g0)

  # Load masks
  masks <- readRDS(masks)

  # Initilize output
  out <- list()

  # Loop through a vector with names of imputed files
  for(i in seq_along(imputed)){

    # Load imputed data for current mask
    g1 <- ReadFastPHASE(imputed[i])

    # Convert to matrix
    c1 <- seq2mat(g1)
    rm(g1)

    # Select mask
    m1 <- masks[[i]]

    # Convert into TRUE/FALSE
    m1 <- m1 == 1

    # Replicate rows
    m1 <- m1[rep(seq_len(nrow(m1)), each = 2),]

    # Subset imputed alleles by mask
    d0 <- c0[m1]
    d1 <- c1[m1]

    # Count statistics
    if(is.null(id)) { out[[i]] <- CountStat(d0, d1)
    } else { out[[i]] <- c(CountStat(d0, d1), id = id) }

  }

  plyr::ldply(out, function(x) x)
}

CountStat <- function(a0, a1){
  # Counts statistics for imputed alleles
  # Args:
  #  a0: original alleles
  #  a1: imputed alleles
  # Returns:
  #  Vector with allele and genotype errors.

  # Check input
  if(length(a0) != length(a1))
    stop(sprintf("Vectors are not equal! Original: %s, imputed: %s",
                 length(a0), length(a1)))

  out <- a0 == a1

  # Count statistics
  nalleles <- length(out)
  ngenotypes <- nalleles/2
  tp <- sum(out) # Number of true positive alleles

  alleles <- (nalleles - tp)/nalleles

  # Count true positive genotypes
  t1 <- which(out == TRUE)
  t2 <- t1[-1]
  dif <- t2 - t1[-tp]

  # How many odd indexes do we have?
  ind <- t1[which(dif == 1)]

  # Count accuracy of genotype imputation
  acc <- sum(is.odd(ind))/ngenotypes

  # Output statistics
  c(alleles = alleles, genotypes = 1 - acc)

}

seq2mat <- function(vec){
  # Converts sequences into matrix

  tmp <- plyr::laply(vec, function(x) strsplit(x, split = ""))
  tmp <- plyr::ldply(tmp, function(x) x)

  tmp <- as.matrix(tmp)
  dimnames(tmp) <- NULL
  tmp

}

# Checks whether the number is odd
is.odd <- function(x) x %% 2 != 0
