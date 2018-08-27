# The sc5mc Package - Analysis of Single Cell Methylation

This package provides a set of tools for the analysis of single cell methylation.
it includes basic pipeline to import data from raw fastq, as well as functions for analysis.

### Code
Source code can be found at: https://bitbucket.org/tanaylab/sc5mc


## Installation
Please make sure that you loaded the MKL module: 


```bash
ml load mkl
```

And then run:


```r
install.packages('sc5mc', repos=c(getOption('repos'), 'https://tanaylab.bitbucket.io/repo'))
```
