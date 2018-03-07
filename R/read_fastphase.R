#'Read fastPHASE files
#'
#'\code{ReadfastPHASE} parses fastPHASE *.inp files. There are two types of
#'them: with ids of samples provided or not. The function loads the second one,
#'simplified. The alleles can be coded as letters or numbers, but missing ones
#'should be as "?".
#'
#'@param inp path/to/filename.inp
#'@return A character vector with strings representing the genotypes. Two
#'  sequencies per individual.
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
  message("Sequences: ", length(l))
  message("Markers: ", nchar(l[1]))

  return(l)
}
