#'Read fastPHASE files
#'
#'Parses fastPHASE *.inp files. There are two types of them: with ids of samples
#'provided or not. The function recognizes the both. The alleles can be coded as
#'letters or numbers, but missing ones should be always "?". Check that alleles
#'in input file are written without whitespaces. The haplotype should look like
#''AACTT' and not like 'A A C T T'.
#'
#'@param inp path/to/filename.inp
#'@return A character vector with strings of letters representing haplotypes.
#'  The length of vector is 2N, where N is the number of diploid individuals.
#'@export
ReadFastPHASE <- function(inp){
  # Check input
if(file.access(inp) == -1) stop(inp, " doesn't exist!", call. = F)

  # Load data set
  message("\nLoading ", inp)
  l <- readLines(inp)

  # Take the number of characters of the last element and
  # leave only the elements that have this number of characters
  l <- l[nchar(l) == nchar(l[length(l)])]

  # Print info about dataset
  message("Haplotypes: ", length(l))
  message("SNPs: ", nchar(l[1]))

  return(l)
}

#' Save a vector with haplotypes as a simplified fastPHASE file ready for
#' imputation.
#'
#' @param g character vector with haplotypes
#' @param pref basename of output file. inp extention is added.
#' @param n ordinal number of masks. If provided, 'm' followed by this number
#'   will be added to the basename of masked files.
#'
#' @return No value
#' @export
WriteFastPHASE <- function(g, pref, n = NULL) {

  # Set output filename
  fn <- ifelse(is.null(n), paste0(pref, ".inp"), sprintf("%s.m%s.inp", pref, n))

  # Remove output file if it exists
  if(file.exists(fn)) file.remove(fn)

  # Get meta data
  N <- length(g) # Number of sequences
  M <- nchar(g[1]) # Number of markers

  # Create output file
  write(c(N/2, M), file = fn, ncolumns = 1)

  # Write masked sequences to file
  plyr::l_ply(g, write, file = fn, append = T)

  message(sprintf("File %s is saved", fn))

}
