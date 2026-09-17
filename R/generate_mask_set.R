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
#' The number of masks in a set can vary in both cases. The masks within one
#' set are drawn disjointly: genotypes hidden by earlier masks are never
#' re-drawn by later ones, so together the masks of a set cover about
#' \code{n * p} of the available genotypes. Call \code{\link{set.seed}}
#' beforehand to make the masks reproducible.
#'
#' @param g Character vector. The elements are sequences made of alleles.
#' The length of \emph{g} equals to 2*\emph{N}, where \emph{N} is the number of
#' individulas, assuming the ploidy of 2.
#' @param n Number of masks to be generated
#' @param p Proportion of genotypes to be masked at each loci or proportion of
#' markers to be masked in chromosome
#' @param type Type of masking. "genotype" for hiding genotypes and "marker" for
#' hiding the markers. The default is "genotype".
#' @param samples Optional character vector of length \emph{N} with sample ids.
#' If given, the ids are attached to the masks as row names.
#' @param markers Optional character vector of length \emph{M} with marker ids
#' (e.g. rs numbers). If given, the ids are attached to the masks as column
#' names. Dimnames make the masks self-describing; they are used by
#' \code{\link{WriteMaskSet}} and enable the consistency checks in
#' \code{\link{EstimateQuality}}.
#'
#' @return A list of length \emph{n} containing masks as binary (0/1) matrices
#'   with \emph{N} rows and \emph{M} columns, where 1 marks a hidden genotype.
#' @export
#'
GenerateMaskSet <- function(g, n, p, type = "genotype", samples = NULL,
                            markers = NULL){

  # Initilize variables
  M <- nchar(g[1]) # number of markers
  N <- length(g)/2 # number of individuals

  if (!is.null(samples) && length(samples) != N)
    stop(sprintf("'samples' has length %s, expected N = %s", length(samples), N))
  if (!is.null(markers) && length(markers) != M)
    stop(sprintf("'markers' has length %s, expected M = %s", length(markers), M))

  dn <- NULL
  if (!is.null(samples) || !is.null(markers))
    dn <- list(if (is.null(samples)) NULL else samples,
               if (is.null(markers)) NULL else markers)

  out <- list()

  if (type == "marker") {
    message("Hiding of markers is choosen")
    # Create first mask
    snps <- 1:M
    size <- floor(p*M)
    ind <- sample(x = snps, size = size)
    out[[1]] <- makeMask(N, M, ind, dn)

    # Make list with other masks
    for (i in 2:n) {
      snps <- snps[!snps %in% ind]
      ind <- sample(x = snps, size = size)
      out[[i]] <- makeMask(N, M, ind, dn)
    }

    return(out)
  }

  message("Hiding of genotypes is choosen")

  # Create a binary matrix with originally missing values
  m0 <- GetMissing(g, M)
  message("Proportion of originally missing genotype: ",
          round(sum(m0)/(nrow(m0) * ncol(m0)), 4))

  # Count available genotypes per marker
  size <- round((N - colSums(m0)) * p)

  # Keep genotypes that are already hidden; updated after each mask so that
  # the masks of one set stay disjoint
  masked <- m0

  # Generate diffferent n masks
  for(i in seq_len(n)) {

    message("Generating mask ", i, "...")

    out[[i]] <- GenerateMask(masked, size, dn)

    # Update genotypes that are already masked
    masked <- masked + out[[i]]
  }
  return(out)
}

GenerateMask <- function(masked, size, dimnames = NULL){
  # Generate a new mask
  # Args:
  #   masked: matrix keeping genotypes already hidden (originally missing or
  #     hidden by the previous masks of the set)
  #   size: vector with the number of genotypes to hide per marker
  #   dimnames: optional dimnames attached to the output mask
  # Returns:
  #   A new mask as a binary matrix

  m <- matrix(0, ncol = length(size), nrow = nrow(masked))

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

  if (!is.null(dimnames)) dimnames(m) <- dimnames
  m
}

GetMissing <- function(data, M){
  # Determines missing genotypes
  # Args:
  #   data: character vector with haplotypes (length 2N)
  #   M: number of markers
  # Returns:
  #   Binary matrix (N x M), where 1 corresponds to missing genotype

  message("Counting missing genotypes...")

  # Split haplotypes into an allele matrix (2N x M)
  x <- matrix(unlist(strsplit(data, "", fixed = TRUE), use.names = FALSE),
              nrow = length(data), byrow = TRUE)

  # A genotype is missing when both of its alleles are unknown
  miss <- (x[seq(1, nrow(x), 2), ] == "?") & (x[seq(2, nrow(x), 2), ] == "?")
  storage.mode(miss) <- "numeric"
  miss
}

makeMask <- function(i, j, x, dimnames = NULL) {
  # Crates matrix with 0, where columns x are filled with 1.
  # Args:
  #  i, j: the number of rows and columns
  #  x: a vector with columns
  # Returns:
  #  Matrix
  m <- matrix(0L, nrow = i, ncol = j)
  m[, x] <- 1L
  if (!is.null(dimnames)) dimnames(m) <- dimnames
  m
}
