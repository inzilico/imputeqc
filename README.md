What is it for
--------------

*imputeqc* is an R package and accompanied scripts to estimate the quality of imputation of genotypes that was made with [fastPHASE](http://scheet.org/software.html) tool. It allows to choose the appropriate number of haplotype clusters (K) for the search of selection fingerprints with [hapFLK](https://forge-dga.jouy.inra.fr/projects/hapflk) test.   

Ho to install from GitHub
--------------------------

1. Make sure you have [devtools](https://github.com/r-lib/devtools) package installed. Run from R

```
install.packages("devtools")
```

2. Install *imputeqc*

```
devtools::install_github("inzilico/imputeqc", build_vignettes = TRUE)
```

How to use
----------

Read a vignette **How to Select the Number of Clusters for fastPHASE**:

* [download *.pdf](vignettes/k_selection.pdf). 
* On a local machine, the vignette can be accessed as follow: 
```
browseVignettes("imputeqc")
```    
* On remote machine, the vignette can be opened in the "Help" tab of RStudio:
```
vignette("k_selection")
```

Contacts
--------
Gennady Khvorykh, a bioinformatician, [inzilico.com](http://inzilico.com)

Interested in contributing to the project? Suggestions, questions, and comments are open! Feel free [to drop me the message](http://www.inzilico.com/contacts/).
