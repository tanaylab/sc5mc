
# The sc5mc Package - Analysis of Single Cell Methylation

This package provides a set of tools for the analysis of single cell
methylation. it includes basic pipeline to import data from raw fastq,
as well as functions for analysis.

### Code

Source code can be found at: <https://github.com/tanaylab/sc5mc>

## Installation

Please make sure that Intel Math Kernel Library (MKL) is installed.

At the tanaylab cluster, you can just run:

``` bash
ml load mkl
```

And then run:

``` r
remotes::install_github("tanaylab/sc5mc")
```
