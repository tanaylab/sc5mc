% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/avg_meth.R
\name{sc5mc.global_meth_trend}
\alias{sc5mc.global_meth_trend}
\title{Global methylation trend per cell}
\usage{
sc5mc.global_meth_trend(
  smat,
  breaks = seq(0, 0.3, by = 0.002),
  include.lowest = TRUE,
  min_cpgs = 1,
  parallel = TRUE
)
}
\arguments{
\item{smat}{smat object}

\item{breaks}{breaks of cpg content}

\item{include.lowest}{if 'TRUE', the lowest value of the range determined by
breaks is included}

\item{min_cpgs}{minimal number of CpGs per cg content strat}

\item{parallel}{compute parallely per group (using doMC package)}
}
\value{
data frame with the following fields:
`cg_cont` with the CpG content bin,
`breaks_numeric` with the middle point between each bin
`cell_id` cell name
`ncpgs` number of CpGs in the bin
`meth` number of methylated CpGs
`unmeth` number of unmethylated CpGs
`avg` fraction of methylated calls out of all calls
}
\description{
Global methylation trend per cell
}
