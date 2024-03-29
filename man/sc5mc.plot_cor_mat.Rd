% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_cells.R
\name{sc5mc.plot_cor_mat}
\alias{sc5mc.plot_cor_mat}
\title{plot cell cell correlation matrix}
\usage{
sc5mc.plot_cor_mat(
  cell_cor,
  row_ord = NULL,
  col_ord = NULL,
  show_colnames = FALSE,
  show_rownames = FALSE,
  breaks = NULL,
  color_pal = NULL,
  ...
)
}
\arguments{
\item{cell_cor}{cell cell correlation (output of sc5mc.calc_pdiff_cor)}

\item{row_ord}{specific order for the rows}

\item{col_ord}{specific order for the columns}

\item{show_colnames}{show column names}

\item{show_rownames}{show row names}

\item{color_pal}{color pallete}

\item{...}{other parameters of pheatmap1}
}
\description{
plot cell cell correlation matrix
}
