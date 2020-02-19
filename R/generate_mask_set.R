#' Generate a set of masks
#'
#' Two types of masks can be generated. The first one hides the genotypes and
#' the second one hides the markers. Actually, both hides the genotypes,
#' but in the first case they are randomly distributed all over the chromosome,
#' while in the second one they are groupped to one or several markers in such a way
#' that the whole marker is hidden. That is why we say that we hide markers.
#'
#' Thus, given a proportion of genotypes to be masked at each loci or
#' a proportion of markers to be masked in chromosome
#' \emph{GenerateMaskSet} samples a set of different masks.
#'
#' The number of masks in a set can vary in both cases.
#'
#' @param g Character vector. The elements are sequences made of alleles.
#' The length of \emph{g} equals to 2*\emph{N}, where \emph{N} is the number of
#' individulas, assuming the ploidy of 2.
#' @param n Number of masks to be generated
#' @param p Proportion of genotypes to be masked at each loci or proportion of
#' markers to be masked in chromosome
#' @param type Type of masking. "genotype" for hiding genotypes and "marker" for
#' hiding the markers. The default is "genotype".
#'
#' @return A list of length \emph{n} containing masks as matrices
#' @export
#'
GenerateMaskSet <- function(g, n, p, type = "genotype"){

  # Initilize variables
  M <- nchar(g[1]) # number of markers
  N <- length(g)/2 # number of individuals
  out <- list()

  if (type == "marker") {
    message("Hiding of markers is choosen")
    # Create first mask
    snps <- 1:M
    size <- floor(p*M)
    ind <- sample(x = snps, size = size)
    out[[1]] <- makeMask(N, M, ind)

    # Make list with other masks
    for (i in 2:n) {
      snps <- snps[!snps %in% ind]
      ind <- sample(x = snps, size = size)
      out[[i]] <- makeMask(N, M, ind)
    }

    return(out)
  }

  message("Hiding of genotypes is choosen")

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
  return(out)
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

makeMask <- function(i, j, x) {
  # Crates matrix with 0, where columns x are filled with 1.
  # Args:
  #  i, j: the number of rows and columns
  #  x: a vector with columns
  # Returns:
  #  Matrix
  m <- matrix(0L, nrow = i, ncol = j)
  m[, x] <- 1
  m
}

