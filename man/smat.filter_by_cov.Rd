% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smat.R
\name{smat.filter_by_cov}
\alias{smat.filter_by_cov}
\title{Filter smat object by cpgs / cells coverage}
\usage{
smat.filter_by_cov(
  smat,
  min_cpgs = 1,
  max_cpgs = Inf,
  min_cells = 1,
  max_cells = Inf
)
}
\arguments{
\item{smat}{smat object}

\item{min_cpgs}{minimal number of CpGs}

\item{max_cpgs}{maximal number of CpGs}

\item{min_cells}{minimal number of cells}

\item{max_cells}{maximal number of cells}
}
\description{
Filter smat object by cpgs / cells coverage
}
