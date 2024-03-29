% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_NNs.R
\name{get_NNs}
\alias{get_NNs}
\title{get_NNs}
\usage{
get_NNs(
  data,
  col_dims,
  method = "FULL",
  R = NA,
  k = NA,
  mode = "hybrid",
  i_omit = c(),
  force_FULL = FALSE,
  multiplier = 1,
  coverage = 1,
  seed = NA,
  verbose = FALSE
)
}
\arguments{
\item{data}{data table containing all necessary columns and rows}

\item{col_dims}{character vector: name of the dimension column(s)}

\item{method}{character vector: method of sampling points to build neighborhoods for. One of:
\itemize{
\item 'FULL' = builds neighborhoods for every points in the dataset (with non-NA data). See \code{\link[=get_NNs.FULL]{get_NNs.FULL()}} for more details.
\item 'FAST' = builds neighborhoods for a much smaller but representative sample of the data's dimensional space. Algorithm can be tweaked further (see below). See \code{\link[=get_NNs.FAST]{get_NNs.FAST()}} for more details.
}}

\item{R}{(optional for mode='k') numeric: radius from each point to build a neighborhood from. Only needed if using 'fr' or 'hybrid' mode}

\item{k}{(optional for mode='fr') integer: number of closest neighbors (including self) to build neighborhood from. Only needed if using 'k' or 'hybrid' mode}

\item{mode}{(optional) character: mode of building neighborhoods. One of:
\itemize{
\item 'fr' = fixed-radius mode. Each anchor's neighborhood is composed of every point within a radius R of the anchor point.
\item 'k' = k-nearest mode. Each anchor's neighborhood is composed of the closest k points to the anchor point.
\item 'hybrid' = fr+k mode. Runs fixed-radius mode and then runs k-nearest mode on neighborhoods with sizes less than k. Recommended for reducing noise.
}}

\item{i_omit}{(optional) numeric vector: indices of data rows for which to not consider at all in nearest neighbor algorithm}

\item{force_FULL}{(optional, method='FULL') logical: whether to continue with function if data's sample size is large (N>20,000). It is recommended to use \code{\link[=locPGSacc.FAST]{locPGSacc.FAST()}} in such situations.}

\item{multiplier}{(optional, method='FAST') integer: number of times each point needs to be assigned to a neighborhood before algorithm ends. Bigger = slower = more anchors.}

\item{coverage}{(optional, method='FAST') numeric: proportion of points that need to meet the multiplier parameter before algorithm ends. Bigger = slower = more anchors}

\item{seed}{(optional, method='FAST') numeric: seed used for algorithm}

\item{verbose}{(optional, method='FAST') logical: whether to print out more detailed algorithm-related information}
}
\value{
Returns a nested list, where the i'th element contains a numeric vector of
the i'th point's neighbors (including element i itself) if computed
}
\description{
Returns nested list of samples' nearest neighbors
}
\details{
This function takes a dataset containing columns for a phenotype,
a PGS for the phenotype, and coordinates for some dimensional-space (e.g.
genetic principal component space) and computes neighborhoods. A point's
neighborhood can be defined in different ways, but broadly consists of points
close to the original point chosen (called the anchor point).
}
\author{
Nuno R. G. Carvalho: \email{nunocarvalho@gatech.edu}
}
