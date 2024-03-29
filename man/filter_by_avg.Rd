% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cgdb.R
\name{filter_by_avg}
\alias{filter_by_avg}
\title{Filter cgdb object by cpgs average methylation}
\usage{
filter_by_avg(db, min_avg = 0, max_avg = 1, min_cov = 1, max_cov = Inf)
}
\arguments{
\item{db}{smat object}

\item{min_avg}{minimal average methylation per CpG}

\item{max_avg}{maximal average methylation per CpG}

\item{min_cov}{minimal coverage per CpG}

\item{max_cov}{maximal coverage per CpG}
}
\description{
Filter cgdb object by cpgs average methylation
}
