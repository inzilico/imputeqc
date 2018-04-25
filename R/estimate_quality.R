#' Estimate the quality of imputation
#'
#' Counts the discordance of genotypes imputed.
#'
#' @param origin path to original unmasked file
#' @param masks path/to/masks.RDS, where \code{masks.RDS} keeps the masks
#'   created upstream with \code{\link{GenerateMaskSet}} and saved as
#'   \code{*.RDS}.
#' @param imputed vector of paths to files imputed. In the case of fastPHASE
#'   program the elements are \code{path/to/*_genotypes.out}. In the case of
#'   BEAGLE - \code{path/to/*.vcf} or \code{path/to/*.vcf.gz}. The order of
#'   elements of vector and masks should coincide.
#' @param id character or numeric, id of computational experiment. In case you
#'   run several calculations with different model parameter to find the best
#'   one, you can mark each run with id for the convinience of further
#'   visualization. The argument is optional.
#'
#' @return A data frame with three columns: "alleles", "genotypes", and "id" if
#'   provided. The first two contains the errors, the third one - id of the
#'   experiment. The error is counted as \code{(1 - accuracy)}, where
#'   \code{accuracy} is a proportion of correctly imputed genotypes and alleles.
#'   By correctly imputed genotype we mean that both alleles coincided with the
#'   original ones. The function returns the values for one set of test files.
#' @export
#'
EstimateQuality <- function(origin, masks, imputed, id = NULL){

  # Load unmasked data set
  h0 <- LoadHaps(origin)

  # Convert to matrix
  c0 <- seq2mat(h0)
  rm(h0)

  # Output missingness for unmasked data
  ms <- c0 == "?"
  message("Missingness of original data: ", round(sum(ms)/length(ms), 4))

  # Load masks
  masks <- readRDS(masks)

  # Initilize output
  out <- vector("list", length(imputed))

  # Loop through all imputed files
  for(i in seq_along(imputed)){

    # Load imputed data as haplotypes
    h <- LoadHaps(imputed[i])

    # Convert haplotypes into matrix with alleles
    c1 <- seq2mat(h)
    rm(h)

    # Select mask
    m <- masks[[i]]

    # Convert into TRUE/FALSE
    m <- m == 1

    # Replicate rows
    m <- m[rep(seq_len(nrow(m)), each = 2),]

    # Subset original and imputed alleles by mask
    a0 <- c0[m]
    a1 <- c1[m]

    # Generate vector with odd and even numbers
    ind1 <- seq(1, length(a0), 2)
    ind2 <- seq(2, length(a0), 2)

    # Make data frame with odd and even numbers
    indx <- data.frame(odd = ind1, even = ind2, stringsAsFactors = F)

    # Count discordance of genotypes
    ct <- vector("logical", nrow(indx))

    for (j in seq_len(nrow(indx))) {
      ct[j] <- sum(a0[c(indx$odd[j], indx$even[j])] %in%
                     a1[c(indx$odd[j], indx$even[j])]) != 2
    }

    disc <- sum(ct)/length(ct)
    disc <- round(disc, 6)
    message("Discordance: ", disc)

    # Make the output
    if(is.null(id)) { out[[i]] <- c(disc = disc)
    } else { out[[i]] <- c(discordance = disc, id = id) }

  }

  plyr::ldply(out, function(x) x)
}

#' Boxpolt discordance
#'
#' Draw boxplots of discordance estimated for different model parameters.
#'
#' @param fname path to dataframe with two columns 'disc' and 'id'
#' @param tl optional title of the plot
#' @param stl optional subtitle of the plot
#'
#' @return ggplot object
#' @import ggplot2
#' @export
#'
PlotDiscordance <- function(fname, tl = NULL, stl = NULL, id = NULL, ...){

  if(!file.exists(fname)) stop((sprintf("File %s doesn't exist", fname)))

  # Load data set
  df <- utils::read.table(fname, header = T)

  ggplot(df, aes(y = discordance, x = id, group = id)) + geom_boxplot(varwidth = TRUE) +
    labs(title = ifelse(is.null(tl),"",tl), subtitle = ifelse(is.null(stl), "", stl),
         y = "Discordance", x = ifelse(is.null(id), "", id)) +
    scale_x_discrete(limits = c(5, 10, 15, 20, 25, 30, 35, 40, 45))

}

seq2mat <- function(vec){
  # Converts sequences into matrix

  tmp <- plyr::laply(vec, function(x) strsplit(x, split = ""))
  tmp <- plyr::ldply(tmp, function(x) x)

  tmp <- as.matrix(tmp)
  dimnames(tmp) <- NULL
  tmp

}

# Checks whether the number is odd
is.odd <- function(x) x %% 2 != 0

LoadHaps <- function(input){
  # Load data and convert them into haplotypes. fastPHASE as well as vcf or vcf
  # gzipped formats are accessable.
  # Args:
  #  input: path/to/filename.{inp,vcf,vcf.gz}
  # Returns: Vector with haplotypes

  # Get extention of input file
  m <- regexpr("\\.([[:alnum:]]+)$", input)
  ext <- regmatches(input, m)

  # Load data
  g <- switch(ext,
              .vcf = GetHaps(input),
              .gz = GetHaps(input),
              .inp = ReadFastPHASE(input),
              .out = ReadFastPHASE(input))

  if(is.null(g)) stop("Haplotypes aren't loaded!", call. = F)

  return(g)

}
