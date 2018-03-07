% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_quality.R
\name{EstimateQuality}
\alias{EstimateQuality}
\title{Estimate the quality of imputation}
\usage{
EstimateQuality(origin, masks, imputed, K)
}
\arguments{
\item{origin}{path/to/filename.inp, where \code{filename.inp} is the
original (unmasked) fastPHASE file}

\item{masks}{path/to/masks.RDS, where \code{masks.RDS} keeps the masks
created upstream with \code{\link{GenerateMaskSet}} and saved as \code{*.RDS}.}

\item{imputed}{vector of \code{path/to/*_genotypes.out}. The files
have genotypes imputed with fastPHASE. They should be in the same
order as masks were generated.}

\item{K}{The number of clusters in the hapFLK test}
}
\value{
A data frame with three columns: "alleles", "genotypes", and "K". The
  first two contains the errors, the third one - \code{K} (the number of
  clusters). The error is counted as \code{1 - accuracy}, where
  \code{accuracy} is a proportion of correctly imputed genotypes and alleles.
  By correctly imputed genotype we mean that both alleles coincided with the
  original ones. The function returns the values for one set of test files.
}
\description{
\code{EstimateQuality} counts the proportion of correctly imputed alleles and
genotypes.
}