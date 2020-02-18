#' Generate a set of masks
#'
#' Given a proportion of genotypes to be masked at each loci
#' \emph{GenerateMaskSet} samples a set of different masks.
#' The number of masks in a set can vary.
#'
#' @param g Character vector. The elements are sequences made of alleles.
#' The length of \emph{g} equals to 2*\emph{N}, where \emph{N} is the number of
#' individulas, assuming the ploidy of 2.
#' @param n Number of masks to be generated
#' @param p Proportion of genotypes to be masked at each loci
#'
#' @return A list of length \emph{n} containing masks as matrices
#' @export
#'
GenerateMaskSet <- function(g, n, p){

  # Initilize variables
  M <- nchar(g[1]) # number of markers
  N <- length(g)/2 # number of individuals
  out <- list()

  # Create a binary matrix with originally missing values
  m0 <- GetMissing(g, M)

  # Count available genotypes per marker
  size <- round((N - colSums(m0)) * p)

  # To keep genotypes that are already masked
  masked <- m0

  # Generate diffferent n masks
  for(i in seq_len(n)) {

    # Initiate empty mask
    m <- matrix(0, ncol = M, nrow = N)

    message(sprintf("Generating mask %s...", i))

    # Generate new mask
    out[[i]] <- GenerateMask(m, masked, size)

    # Update genotypes that are alreade masked
    masked <- Reduce('+', out)
  }
  out
}

GenerateMask <- function(m, masked, size){
  # Generate a new mask
  # Args:
  #   m: matrix with zeros
  #   masked: matrix keeping genotypes previously masked
  #   size: vector with the number of genotypes available for masking
  # Returns:
  #   A new mask as a binary matrix

  for (j in seq_len(ncol(m))) {

    ind <- which(masked[, j] == 0)

    # If no genotypes left for masking skip this position
    if(length(ind) == 0 || size[j] == 0) next

    # Sample
    if(length(ind) >= size[j]) {
      add <- sample(ind, size[j])
    } else { add <- ind }

    m[add, j] <- 1
  }
  m
}


GetMissingMarkers <- function(data){
  # Determines the positions of missing bases per individual
  # Args:
  #   data: character vector of sequnces
  # Returns:
  #   List with positions of missed genotypes for each individual
  mindex <- sapply(data, function(x) unlist(gregexpr("\\?", x)), USE.NAMES = F)
  mindex[seq(2, length(mindex), 2)]
}

ToBinary <- function(data, jmax){
  # Converts list with the positions into binary matrix
  # Args:
  #  data: list with the positions
  #  jmax: number of markers
  # Returns:
  #  Binary matrix, where 1 means missing and 0 - nonmissing values.
  sapply(seq_len(jmax), function(j){
    sapply(data, function(i) ifelse(j %in% i, 1, 0))
  })
}

GetMissing <- function(data, M){
  # Determines missing genotypes
  # Args:
  #   data: character vector with sequences
  #   M: number of markers
  # Returns:
  #   Binary matrix, where 1 corresponds to missing genotype.

  message("Counting missing genotypes...")

  # Get list of missing markers per individual
  tmp <- GetMissingMarkers(data)

  # Convert list of missing markers into binary matrix
  m <- ToBinary(tmp, M)

  # Print proportion of missing genotypes
  p <- sum(colSums(m))/(dim(m)[1] * dim(m)[2])
  message("Proportion of originally missing genotype: ", round(p, 4))

  m
}


