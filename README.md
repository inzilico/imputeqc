What `imputeqc` is for
----------------------

`imputeqc` is an R package and accompanied scripts to estimate the quality of imputation of genotypes that was made with [fastPHASE](http://scheet.org/software.html) tool. It allows to choose the appropriate number of haplotype clusters (K) for the search of selection fingeprints with the use of [hapFLK](https://forge-dga.jouy.inra.fr/projects/hapflk) test.   

Ho to install a package from GitHub
-----------------------------------

1. Make sure you have [devtools](https://github.com/r-lib/devtools) package installed. Run from R

    install.packages("devtools")
  
2. Load the `devtools` package and install `imputeqc` from GitHub

    library("devtools")
    install.github("inzilico/imputeqc")

How to use a package
--------------------

Check a vignette from the package. Run in R 

    browseVignettes("imputeqc")

Contacts
--------
Gennady Khvorykh, a bioinformatician, [inzilico.com](http://inzilico.com)

Interested in contributing to the project? Suggestions, questions, and comments are open! Feel free [to drop me the message](http://www.inzilico.com/contacts/).
