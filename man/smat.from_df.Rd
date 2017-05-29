% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smat.R
\name{smat.from_df}
\alias{smat.from_df}
\title{Create smat object from data frame}
\usage{
smat.from_df(df)
}
\arguments{
\item{df}{intervals (chrom, start, end fields) data frame with additonal fields:
meth (methylated calls), unmeth (unmethylated calls) and cov (total coverage)}
}
\value{
smat object
}
\description{
Create smat object from data frame
}