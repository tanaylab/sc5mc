% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cgdb.R
\name{count_pairs}
\alias{count_pairs}
\title{Count methylation calls (00,01,10,11) for pairs of cells}
\usage{
count_pairs(db, max_chunk_size = 200)
}
\arguments{
\item{db}{cgdb object}
}
\value{
data frame with the following fields: cell1, cell2, n00, n01, n10, n11
}
\description{
Count methylation calls (00,01,10,11) for pairs of cells
}
