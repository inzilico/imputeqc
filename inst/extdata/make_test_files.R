#!/usr/bin/Rscript
# Script reads inp (fastPHASE input file) or vcf, generates masks to hide randomly a
# given proportion of genotypes or a given proportion of markers, applies masks to the input file, and saves test
# files thus produced. Test files can be further sent to fastPHASE and BEAGLE for
# imputation.
# Author: Gennady Khvorykh, http://inZilico.com

library(imputeqc)
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
                metavar = "outputpath"),
    make_option(c("-t", "--type"), type = "character", default = "genotype",
                help = "Type of masking (genotype/marker). To hide genotypes set 'genotype',
                and to hide marker set 'marker'. [default %default]",
                metavar = "type"),
    make_option(c("-s", "--swap"), type = "character", default = "F",
                help = "If T, the haplotypes of diploid organism are swapped [default %default]",
                metavar = "swap")
  )

  # Create parser
  parser <- OptionParser(usage = "%prog [options] input",
                         option_list = option.list,
                         prog = "make_test_files.R",
                         description = "\ninput: full/path/to/filename.{inp,vcf}")

  # Read arguments from command line
  arguments <- parse_args(parser, positional_arguments = TRUE)

  # Check input file
  input <- arguments$args

  if (length(input) != 1) {
    cat("\nIncorrect number of required positional arguments!\n\n")
    return(print_help(parser))
  }

  # Check input file exist
  if(!file.exists(input)) stop(input, " doesn't exist!", call. = F)

  # Get option arguments
  options <- arguments$options
  output <- options$output
  swap <- options$swap

  # To measure the computation time
  start <- Sys.time()

  # Get extention of input file
  m <- regexpr("\\.([[:alnum:]]+)$", input)
  ext <- regmatches(input, m)

  # Load data
  if(ext == ".vcf" | ext == ".gz") {
    data <- ReadVCF(input, swap)
    g <- data[["haps"]]
    vcf <- data[["vcf"]]
  }

  if(ext == ".inp") {
    g <- ReadFastPHASE(input)
    vcf <- NULL
  }

  # Generate n masks
  if(is.null(g)) stop("Haplotypes aren't loaded. Probably unknown type of input file!",
                      call. = F)

  masks <- GenerateMaskSet(g = g,
                           n = options$ntest,
                           p = options$proportion,
                           type = options$type)

  # Save `masks` as `masks.RDS` at the same directory as `output`
  dir <- dirname(output)
  if(!dir.exists(dir)) dir.create(dir)
  fn <- paste0(dir, "/masks.RDS")
  saveRDS(masks, fn)
  message(sprintf("File %s is saved", fn))

  # Apply masks to sequences and save the result in files
  ApplyMasks(g, masks, output, vcf)

  # Show time elapsed
  Sys.time() - start

} # End of main()

main()
