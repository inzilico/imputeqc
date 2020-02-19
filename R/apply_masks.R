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
#' @importFrom plyr progress_text
#'
ApplyMasks <- function(g, masks, pref, vcf = NULL) {

  # Initilize
  N <- length(g) # The number of haplotypes

  # Loop through all masks and save the masked data obtained
  for (n in seq_along(masks)) {

    # Get indexes of masked genotypes
    ind <- plyr::alply(masks[[n]], 1, function(v) which(v == 1))

    # Replicate indexes
    ind <- rep(ind, each = 2)

    message(sprintf("Applying mask %s...", n))

    # Loop throug all sequences and mask them
    gm <- plyr::llply(seq_len(N), .progress = plyr::create_progress_bar(name = "text"),
                function(i, g, ind) MaskSequence(g[i], ind[[i]]), g = g, ind = ind)

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

  # Loop throug all positions to be replaced
  for(x in positions) substr(sequence, start = x, stop = x) <- symbol
  sequence
}
