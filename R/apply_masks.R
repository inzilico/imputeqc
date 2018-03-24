#' Apply set of masks to sequencies.
#'
#' Function applies a set of masks to strings obtained with readFastPhase. The output is saved as
#' fastPHASE simpilfied *.inp files. Simplified means that it doesn't include
#' sample ids. The files created are ready for imputation with fastPHASE.
#'
#' @param g Character vector with sequences
#' @param masks List of masks as binary matrices
#' @param pref path/to/prefix to save the result. The number of files created
#'   equals to the length of masks. The filenames are generated
#'   automatically like this: \code{prefix.m{n}.inp}, where \code{prefix} is a user defined
#'   string, \code{n} is an ordinal number of the mask.
#'
#' @return No values
#' @export
#' @importFrom plyr progress_text
#'
#' @examples
ApplyMasks <- function(g, masks, pref) {

  # Initilize
  N <- length(g) # Number of sequences
  M <- nchar(g[1]) # Number of markers

  # Loop throug all masks
  for (n in seq_along(masks)) {

    # Get indexes of masked genotypes
    ind <- plyr::alply(masks[[n]], 1, function(v) which(v == 1))

    # Replicate indexes
    ind <- rep(ind, each = 2)

    message(sprintf("Applying mask %s...", n))

    # Lopp throug all sequences and mask them
    gm <- plyr::llply(seq_len(N), .progress = plyr::create_progress_bar(name = "text"),
                function(i, g, ind) MaskSequence(g[i], ind[[i]]), g = g, ind = ind)

    # Set output filename
    fn <- sprintf("%s.m%s.inp", pref, n)

    # Remove output file if it exists
    if(file.exists(fn)) file.remove(fn)

    # Create output file
    write(c(N/2, M), file = fn, ncolumns = 1)

    # Write masked sequences to file
    plyr::l_ply(gm, write, file = fn, append = T)

    message(sprintf("File %s is saved", fn))
  }

}

MaskSequence <- function(sequence, positions, symbol = "?"){
  # Replace bases in a sequence by `mask` symbol
  # Input:
  #   sequence: character string
  #   positions: vector with positions of characters to be replaced by symbol
  #   symbol: character representing a missing value
  # Returns:
  #   Sequence where some bases are replaced by symbol

  # Loop throug all positions to be replaced
  for(x in positions) substr(sequence, start = x, stop = x) <- symbol
  sequence
}
