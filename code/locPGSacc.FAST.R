# FAST version of locPGSacc
# returns inputted data table with 'n_neighbors' and 'locPGSacc' columns

library(tidyverse)
library(dbscan)

locPGSacc.FAST <- function (
    data,
    col_dims,
    col_pheno,
    col_PGS,
    R = -1,
    k = -1,
    mode = "hybrid", # fr=fixed-radius, k=k closest neighbors, hybrid=highest n_neighbors of fr and k
    multiplier = 1,
    coverage = 1, #between 0 and 1
    seed = NA,
    verbose=FALSE
) {
  # sets seed if given
  if (!is.na(seed)) { set.seed(seed) }
  
  # Checks if col_dims are in data table
  if ( !all(col_dims %in% colnames(data))) {stop("Data table does not contain all specified col_dims")}
  # Checks if data table already contains an output column
  cols_needed <- c("n_neighbors","locPGSacc")
  if ( any(cols_needed %in% colnames(data)) ) {
    cols_duplicated <- cols_needed[cols_needed %in% colnames(data)]
    warning(paste0(c("Data table already contains at least one column used for output purposes:", cols_duplicated), collapse=" "))
  }
  
  # gets number of samples
  data_dims <- data %>% select(any_of(col_dims))
  n_samples <- nrow(data_dims)
  if (k >= n_samples) {warning("k=",k," is bigger than number of samples")}
  
  # empty vectors for later loop
  anchors <- c()
  NN_ids <- list()
  housing <- rep(0,n_samples) # vector of how many housing the i'th element belongs to
  
  if (verbose) {print(paste0("Looping until less than ",ceiling(n_samples*(1-coverage))," unhoused"))}
  
  print(paste0("Computing nearest neighbors for ",n_samples," samples using ",mode,
               "-mode with the settings", ifelse(mode!="k"," R = ",R), ifelse(mode!="fr"," k = ",k)))
  while(sum(housing < multiplier) > n_samples*(1-coverage)) {
    unhoused <- (1:n_samples)[housing < multiplier]
    if (length(unhoused) == 1) {anchor <- unhoused
    } else {anchor <- sample(c(unhoused), 1)}
    
    anchors <- c(anchors, anchor)
    if (verbose) {
      print(paste("Anchor #",length(anchors),"Unhoused:",length(unhoused)))
    }
    
    if (mode=="fr") {
      # Fixed-radius mode
      NNfr <- dbscan::frNN(data_dims, eps=R, sort=FALSE, query=data_dims[anchor,])
      NN_anchor <- NNfr$id[[1]]
    } else if (mode=="k") {
      # k-nearest mode
      NNk <- dbscan::kNN(data_dims, k=k, sort=FALSE, query=data_dims[anchor,])
      NN_anchor <- NNk$id[1,] %>% unname()
    } else if (mode=="hybrid") {
      # hybrid fr+k mode
      # computes NN with fixed-radius
      NNfr <- dbscan::frNN(data_dims, eps=R, sort=FALSE, query=data_dims[anchor,])
      NN_anchor <- NNfr$id[[1]]
      # computes NN with k-nearest for points with <k fixed-radius neighbors
      if (length(NN_anchor) < k) {
        NNk <- dbscan::kNN(data_dims, k=k, sort=FALSE, query=data_dims[anchor,])
        NN_anchor <- split(NNk$id,seq_len(nrow(NNk$id)))[[1]]
      }
    }
    housing[ NN_anchor ] <- housing[ NN_anchor ] + 1
    NN_ids[[anchor]] <- NN_anchor
  }
  
  # gets number of neighbors
  data$n_neighbors <- as.numeric(NA)
  data$n_neighbors[anchors] <- sapply(NN_ids, length)[anchors]
  
  # computes correlation between phenotype and PGS for phenotype
  print("Computing correlation within each neighborhood")
  r_values <- c()
  for (i in anchors) {
    r <- cor(data[ NN_ids[[i]], ][[ col_pheno ]],
             data[ NN_ids[[i]], ][[ col_PGS ]])
    r_values[i] <- r
  }
  data$locPGSacc <- as.numeric(NA)
  data$locPGSacc[anchors] <- r_values[anchors]
  
  return(data)
}
