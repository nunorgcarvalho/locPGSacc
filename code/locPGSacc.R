library(tidyverse)
library(data.table)
library(dbscan)

locPGSacc <- function(
    data,
    col_dims,
    col_pheno,
    col_PGS,
    R = -1,
    k = -1,
    fixed_radius=TRUE
    ) {
  
  # Checks if data table already contains an output tha
  cols_needed <- c("n_neighbors","locPGSacc")
  if (any(cols_needed %in% colnames(data))) {
    cols_duplicated <- cols_needed[cols_needed %in% colnames(data)]
    print(cat("ERROR: Data table already contains at least one column used for output purposes:", cols_duplicated))
    stop()
  }
  
  data_dims <- data %>% select(any_of(col_dims))
  n_samples <- nrow(data_dims)
  
  t1 <- Sys.time()
  if (fixed_radius) {
    print(paste0("Computing fixed-radius=",r," nearest neighbors for ",n_samples," samples"))
    NN <- dbscan::frNN(data_dims, eps=R, sort=FALSE)
    NN_ids <- NN$id
    data$n_neighbors <- sapply(NN$id, length)
  } else {
    print(paste0("Computing k=",k," nearest neighbors for ",n_samples," samples"))
    NN <- dbscan::kNN(data_dims, k=k, sort=FALSE)
    NN_ids <- split(NNk$id,seq_len(nrow(NNk$id)))
    data$n_neighbors <- k
  }
  Sys.time() - t1
  
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