What is it
----------

*imputeqc* is an R package and accompanied scripts to estimate the quality of imputation of genotypes that was made with [fastPHASE](http://scheet.org/software.html) and [BEAGLE](https://faculty.washington.edu/browning/beagle/beagle.html) softwares. Any other tools that support *.inp fastPHASE or VCF files can be also considered. The package is based on masked data analysis. 

Possible applications
---------------------

1. Estimation of the error of gynotype imputation.

2. Optimization of the imputation model parameters, e.g., the number of haplotype clusters. The parameter can be further used for the search of signatures of selection with [hapFLK](https://forge-dga.jouy.inra.fr/projects/hapflk) test.

3. Testing different reference panels for imputation.

4. Benchmarking of different imputation softwares and strategies.

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

Read a vignette [How to Select the Number of Clusters for fastPHASE](https://htmlpreview.github.io/?https://github.com/inzilico/imputeqc/blob/master/vignettes/k_selection.html). 

* On a local machine, the vignette can be accessed as follow: 
```
browseVignettes("imputeqc")
```    
* On remote machine, the vignette can be opened in the "Help" tab of RStudio:
```
vignette("k_selection")
```

License
-------
[MIT](https://en.wikipedia.org/wiki/MIT_License)

Citing
------
Gennady Khvorykh. inzilico/imputeqc v1.0.0 (2018). GitHub repository, https://github.com/inzilico/imputeqc.

Contacts
--------
Gennady Khvorykh, a bioinformatician, [inzilico.com](http://inzilico.com)

Interested in contributing to the project? Suggestions, questions, and comments are open! Feel free [to drop me the message](http://www.inzilico.com/contacts/).
