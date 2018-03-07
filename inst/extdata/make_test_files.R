#!/usr/bin/Rscript
# Script reads fastPHASE input file (*.inp), generates masks to hide randomly a
# given proportion of genotypes, applies masks to the input file, and saves test
# files thus produced. Test files can be further sent to fastPHASE for
# imputation.
# Author: Gennady Khvorykh, http://inZilico.com

suppressPackageStartupMessages(library("optparse"))

main <- function(){

  # Set option arguments
  option.list <- list(
    make_option(c("-n", "--ntest"), type = "integer", default = as.integer(5),
                help = "Number of test files to be generated [default %default]",
                metavar = "n"),
    make_option(c("-p", "--proportion"), type = "double", default = 0.1,
                help = "Proportion of missed genotypes [default %default]",
                metavar = "p"),
    make_option(c("-o", "--output"), type = "character", default = "test/test",
                help = "output/path/prefix to save test files [default %default]",
                metavar = "outputpath")
  )

  # Create parser
  parser <- OptionParser(usage = "%prog + [options] input",
                         option_list = option.list,
                         prog = "make_test_files.R",
                         description = "\ninput: full/path/to/filename.inp, where \"filename.inp\" is fastPHASE input file")

  # Read arguments from command line
  arguments <- parse_args(parser, positional_arguments = TRUE)

  # Check input file
  input <- arguments$args

  if (length(input) != 1) {
    cat("\nIncorrect number of required positional arguments!\n\n")
    return(print_help(parser))
  }

  # Get option arguments
  options <- arguments$options
  output <- options$output

  # To measure the computation time
  start <- Sys.time()

  # Read fastPHASE input file
  g <- ReadFastPHASE(input)

  # Generate n masks
  masks <- GenerateMaskSet(g, n = options$ntest, p = options$proportion)

  # Save `masks` as `masks.RDS` at the same directory as `output`
  dir <- dirname(output)
  if(!dir.exists(dir)) dir.create(dir)
  fn <- paste0(dir, "/masks.RDS")
  saveRDS(masks, fn)
  message(sprintf("File %s is saved", fn))

  # Apply masks to sequences and save the result in files
  ApplyMasks(g, masks, output)

  # Show time elapsed
  Sys.time() - start

} # End of main()

main()
