library(tidyverse)
library(data.table)
library(dbscan)

locPGSacc <- function (
    data,
    col_dims,
    col_pheno,
    col_PGS,
    R = -1,
    k = -1,
    mode = "fr" # fr=fixed-radius, k=k closest neighbors, hybrid=highest n_neighbors of fr and k
    ) {
  
  # Checks if data table already contains an output tha
  cols_needed <- c("n_neighbors","locPGSacc")
  if (any(cols_needed %in% colnames(data))) {
    cols_duplicated <- cols_needed[cols_needed %in% colnames(data)]
    stop(paste0(c("Data table already contains at least one column used for output purposes:", cols_duplicated), collapse=" "))
  }
  
  data_dims <- data %>% select(any_of(col_dims))
  n_samples <- nrow(data_dims)
  if (k >= n_samples) {warning("k=",k," is bigger than number of samples")}
  
  t1 <- Sys.time()
  if (mode=="fr") {
    # Fixed-radius mode
    print(paste0("Computing fixed-radius=",R," nearest neighbors for ",n_samples," samples"))
    NN <- dbscan::frNN(data_dims, eps=R, sort=FALSE)
    NN_ids <- NN$id
    data$n_neighbors <- sapply(NN$id, length)
  } else if (mode=="k") {
    # k-nearest mode
    print(paste0("Computing k=",k," nearest neighbors for ",n_samples," samples"))
    NN <- dbscan::kNN(data_dims, k=k, sort=FALSE)
    NN_ids <- split(NNk$id,seq_len(nrow(NNk$id)))
    data$n_neighbors <- k
  } else if (mode=="hybrid") {
    # hybrid fr+k mode
    print(paste0("Computing hybrid fixed-radius=",R," & k=",k," nearest neighbors for ",n_samples," samples"))
    # computes NN with fixed-radius
    NNfr <- dbscan::frNN(data_dims, eps=R, sort=FALSE)
    NN_ids <- NNfr$id
    # computes NN with k-nearest for points with <k fixed-radius neighbors
    small_neighborhoods <- (1:nrow(data))[sapply(NNfr$id, length) < k]
    NNk <- dbscan::kNN(data_dims, k=k, sort=FALSE, query=data_dims[small_neighborhoods,])
    NN_ids[small_neighborhoods] <- split(NNk$id,seq_len(nrow(NNk$id)))
    data$n_neighbors <- sapply(NN_ids, length)
  }
  Sys.time() - t1
  
  print("Computing correlation within each neighborhood")
  r_values <- c()
  # computes correlation between phenotype and PGS for phenotype
  for (i in 1:nrow(data)) {
    r <- cor(data[ NN_ids[[i]], ][[ col_pheno ]],
             data[ NN_ids[[i]], ][[ col_PGS ]])
    r_values[i] <- r
  }
  data$locPGSacc <- r_values
  
  return(data)
}