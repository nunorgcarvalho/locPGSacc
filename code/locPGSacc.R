# returns inputted data table with 'n_neighbors' and 'locPGSacc' columns

library(tidyverse)
library(dbscan)
source("add_self2neighborhood.R")

locPGSacc <- function (
    data,
    col_dims,
    col_pheno,
    col_PGS,
    R = -1,
    k = -1,
    mode = "hybrid", # fr=fixed-radius, k=k closest neighbors, hybrid=highest n_neighbors of fr and k
    force_largeN = FALSE
) {
  
  # Checks if correct parameter was supplied
  if (mode != "k" & R=-1) {stop("R must be supplied")}
  if (mode != "fr" & k=-1) {stop("k must be supplied")}
  
  # Checks if col_dims are in data table
  if ( !all(col_dims %in% colnames(data))) {stop("Data table does not contain all specified col_dims")}
  # Checks if data table already contains an output column
  cols_needed <- c("n_neighbors","locPGSacc")
  if ( any(cols_needed %in% colnames(data)) ) {
    cols_duplicated <- cols_needed[cols_needed %in% colnames(data)]
    warning(paste0(c("Data table already contains at least one column used for output purposes:", cols_duplicated), collapse=" "))
  }
  
  data_dims <- data %>% select(any_of(col_dims))
  n_samples <- nrow(data_dims)
  if (n_samples > 20000 & !force_largeN) {
    stop("N is (probably) too big for this algorithm. Use locPGSacc.FAST() instead,
         or set 'force_largeN'to TRUE if you are sure you want to proceed")
    }
  if (k >= n_samples) {warning("k=",k," is bigger than number of samples")}
  
  print(paste0("Computing nearest neighbors for ",n_samples," samples using ",mode,
               "-mode with the settings",
               ifelse(mode!="k",paste(" R =",R),""),
               ifelse(mode!="fr",paste(" k =",k),"")))
  if (mode=="fr") {
    # Fixed-radius mode
    NNfr <- dbscan::frNN(data_dims, eps=R, sort=FALSE)
    NN_ids <- NNfr$id %>% add_self2neighborhood()
    data$n_neighbors <- sapply(NN_ids, length)
  } else if (mode=="k") {
    # k-nearest mode
    NNk <- dbscan::kNN(data_dims, k=k, sort=FALSE)
    NN_ids <- split(NNk$id,seq_len(nrow(NNk$id))) %>% add_self2neighborhood()
    data$n_neighbors <- k
  } else if (mode=="hybrid") {
    # hybrid fr+k mode
    # computes NN with fixed-radius
    NNfr <- dbscan::frNN(data_dims, eps=R, sort=FALSE)
    NN_ids <- NNfr$id
    # computes NN with k-nearest for points with <k fixed-radius neighbors
    small_neighborhoods <- (1:nrow(data))[sapply(NNfr$id, length) < k]
    if (length(small_neighborhoods > 0)) {
      NNk <- dbscan::kNN(data_dims, k=k, sort=FALSE, query=data_dims[small_neighborhoods,])
      NN_ids[small_neighborhoods] <- split(NNk$id,seq_len(nrow(NNk$id)))
    }
    NN_ids <- add_self2neighborhood(NN_ids)
    data$n_neighbors <- sapply(NN_ids, length)
  }
  
  # computes correlation between phenotype and PGS for phenotype
  print("Computing correlation within each neighborhood")
  r_values <- c()
  for (i in 1:nrow(data)) {
    r <- cor(data[ NN_ids[[i]], ][[ col_pheno ]],
             data[ NN_ids[[i]], ][[ col_PGS ]])
    r_values[i] <- r
  }
  data$locPGSacc <- r_values
  
  return(data)
}