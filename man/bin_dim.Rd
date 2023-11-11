% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_bandPGS_decay.R
\name{bin_dim}
\alias{bin_dim}
\title{bin_dim}
\usage{
bin_dim(data, ref_window = 0.95, bands = 15)
}
\arguments{
\item{ref_window}{(optional) numeric: proportion of the distribution, centered at the median, to include when making bands}

\item{bands}{(optional) integer: number of bands (i.e. sections) to split up the middle of the sample by according to the dimensional variable. An extra band is made on each end to include dispersed samples.}
}
\value{
Returns a vector of dim groups to apply directly to the dataset:
}
\description{
splits a complete, renamed dataset by bins along some variable
}
\details{
Takes a dataset with complete cases and with column names standardized
('dim','pheno','PGS') and splits individuals into bins depending on their 'dim'
variable. Internal function.
}
\author{
Nuno R. G. Carvalho: \email{nunocarvalho@gatech.edu}
}