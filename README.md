The sc5mc Package - Analysis of Single Cell Methylation
=======================================================

This package provides a set of tools for the analysis of single cell methylation. it includes basic pipeline to import data from raw fastq, as well as functions for analysis.

### Code

Source code can be found at: <https://bitbucket.org/tanaylab/sc5mc>

### Installation

#### Requirements

-   *gpatterns* (<https://bitbucket.org/tanaylab/gpatterns>)

#### Installing gpatterns package:

Download and install *gpatterns*:

``` r
devtools::install_bitbucket("tanaylab/gpatterns", ref='default')
devtools::install_bitbucket("tanaylab/sc5mc", ref='default')
library(sc5mc)
```

#### Using the package

Please refer to the package vignettes for usage and workflow.