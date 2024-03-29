% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smat.R
\name{smat.cell_pairs_marginals}
\alias{smat.cell_pairs_marginals}
\title{Calculate the joint coverage of all pairs of cells}
\usage{
smat.cell_pairs_marginals(smat, cols = NULL)
}
\arguments{
\item{smat}{smat object}

\item{cols}{column names of columns to use. If NULL all columns (cells) would be returned}
}
\value{
tibble with cell pairs (cell1, cell2 fields) and 'ntot' with the
number of jointly covered CpGs
}
\description{
Calculate the joint coverage of all pairs of cells
}
