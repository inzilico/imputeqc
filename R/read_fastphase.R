#'Read fastPHASE files
#'
#'Parses fastPHASE *.inp files. There are two types of them: with ids of samples
#'provided or not. The function loads both types. The alleles can be coded as
#'letters or numbers, but missing ones should be "?".
#'
#'@param inp path/to/filename.inp
#'@return A character vector with strings of letters representing haplotypes.
#'  The length of vector is 2N, where N is the number of diploid individuals.
#'@export
ReadFastPHASE <- function(inp){
  # Check input
if(file.access(inp) == -1) stop(inp, " doesn't exist!", call. = F)

  # Load data set
  message("Loading...")
  l <- readLines(inp)

  # Take the number of characters of the last element and
  # leave only the elements that have this number of characters
  l <- l[nchar(l) == nchar(l[length(l)])]

  # Print info about dataset
  message("Haplotypes: ", length(l))
  message("SNPs: ", nchar(l[1]))

  return(l)
}
