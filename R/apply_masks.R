#' Apply set of masks to haplotypes.
#'
#' Applies a set of masks to haplotypes obtained with readFastPhase. The output
#' is saved as fastPHASE *.inp file not including sample ids (a simplified
#' version of fastPHASE files). It is ready for imputation with fastPHASE.
#'
#' @param g Character vector with original (unmasked) haplotypes
#' @param masks List of masks as binary matrices
#' @param pref path/to/prefix to save the output. The number of files created
#'   equals to the length of masks variable. The filenames are generated
#'   automatically like this: \emph{prefix.m{n}.inp}, where \emph{prefix} is a
#'   user defined string, \emph{n} is an ordinal number of the mask.
#' @param vcf VCF-class object. If provided, the output will be saved as vcf
#'   file. If not, as fastPHASE inp file (default).
#'
#' @return No values
#' @export
#'
ApplyMasks <- function(g, masks, pref, vcf = NULL) {

  # Initilize
  N <- length(g) # The number of haplotypes

  # Loop through all masks and save the masked data obtained
  for (n in seq_along(masks)) {

    # Get indexes of masked genotypes per individual and replicate them
    # for the two haplotypes
    ind <- lapply(seq_len(nrow(masks[[n]])),
                  function(i) which(masks[[n]][i, ] == 1))
    ind <- rep(ind, each = 2)

    message("Applying mask ", n, "...")

    # Loop throug all sequences and mask them
    gm <- vapply(seq_len(N), function(i) MaskSequence(g[i], ind[[i]]),
                 character(1), USE.NAMES = FALSE)

    # Save output
    if(is.null(vcf)) { WriteFastPHASE(gm, pref, n)
      } else { UpdateVCF(gm, pref, n, vcf) }

  }

}

MaskSequence <- function(sequence, positions, symbol = "?"){
  # Replace bases in a sequence by 'symbol' variable
  # Input:
  #   sequence: character string
  #   positions: vector with positions of characters to be replaced by symbol
  #   symbol: character representing a missing value
  # Returns:
  #   Sequence where some bases are replaced by symbol

  if (length(positions) == 0) return(sequence)

  ch <- strsplit(sequence, "", fixed = TRUE)[[1]]
  ch[positions] <- symbol
  paste(ch, collapse = "")
}
