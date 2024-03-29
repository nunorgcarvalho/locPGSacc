% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotPGSdecay.R
\name{plotPGSdecay}
\alias{plotPGSdecay}
\title{plotPGSdecay}
\usage{
plotPGSdecay(
  data,
  col_dim,
  col_pheno,
  col_PGS,
  i_omit = c(),
  ref_window = 0.95,
  bins = 15,
  min_samples = 30,
  fixed_ymin = TRUE,
  show_stats = TRUE,
  plot_rel_dist = FALSE
)
}
\arguments{
\item{data}{data table containing all necessary columns and rows}

\item{col_dim}{character: column name of the dimensional variable for which to compute bins from}

\item{col_pheno}{character: column name of the phenotype of interest}

\item{col_PGS}{character: column name of the polygenic scores for the phenotype of interest}

\item{i_omit}{(optional) numeric vector: indices of data rows for which to not consider at all in nearest neighbor algorithm}

\item{bins}{(optional) integer: number of bins to split up the middle of the sample by according to the dimensional variable. The leftmost and rightmost bins will each contain N/bins individuals. The remaining bins will be of varying sample size but equal range}

\item{min_samples}{(optional) integer: minimum number of samples needed in a bin to compute PGS accuracy. Must be > 2}

\item{fixed_ymin}{(optional) logical: whether the y-axis should start at 0 (default) or just be the minimum accuracy}

\item{show_stats}{(optional) logical: whether the decay slope, standardized decay slope, and corresponding p-value are shown on the plot (top right)}

\item{plot_rel_dist}{(optional) logical: whether the x-axis of the plot should be scaled such that the maximum distance is 1. Does not affect portability statistics.}
}
\value{
Plot of the results of \code{\link[=PGSdecay]{PGSdecay()}}.
}
\description{
Computes PGS accuracy along bins of a continuous variable
}
\details{
Plots the decay of PGS accuracy along some dimensional variables by
cutting up the dimensional variable into bins. More details in \code{\link[=PGSdecay]{PGSdecay()}}.
}
\author{
Nuno R. G. Carvalho: \email{nunocarvalho@gatech.edu}
}
