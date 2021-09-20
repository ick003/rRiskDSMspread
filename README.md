
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rRiskDSMspread

<!-- badges: start -->

[![lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)  
<!-- badges: end -->

Suite of code, and documentation to reproduce the results presented in
the article by Ickowicz et al (2021) modelling the spread and dispersal
of dominant sterile male mosquitoes.

The package is organized as follows:

-   the R directory stores the functions used to run the analysis
-   the inst/code\_article directory stores the different scripts used
    in the analysis. By running the script you reproduce the analysis
-   the inst/extdata directory stores different external datasets used
    in the analysis

More details are provided below.

## Installation

You can install from [GitHub](https://github.com/ick003/rRiskDSMspread)
with:

``` r
## install.packages("remotes") ## if needed
remotes::install_github("ick003/rRiskDSMspread")
```

You also need to install the some utilities package called spreadR (not
to confuse with the CRAN available spreadr) which we make available
under /additional\_resources.

``` r
spreadR_path <- system.file("additional_resources", "spreadR_1.0.27_R_x86_64-pc-linux-gnu.tar.gz", package="rRiskDSMspread")
install.packages(spreadR_path, repos = NULL, type="source")
```

Also, other packages are required for data manipulation and
vizualization:

``` r
list.of.packages <- c(
    "ggplot2", "coda", "abind", "dplyr", "ggmcmc", "scales", "stringr"
  )
  
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
```

## Run the scripts

Once you have installed the package, you can run the scripts (located in
the /code\_article directory) to reproduce the analysis performed in the
article. The scripts are number from 1 to 8:

For example, the prior distributions for the PDE parameters are
constructed using the elicited values from the experts. Check script 01:

``` r
source("/code_article/01_elicited_priors.R")
```

or open the file, and run line by line. And so on with all the scripts.
