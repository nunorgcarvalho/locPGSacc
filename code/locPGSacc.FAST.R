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
    coverage = 1,
    seed = NA,
    verbose=FALSE
) {
  if (!is.na(seed)) { set.seed(seed) }
  
  data_dims <- data %>% select(any_of(col_dims))
  n_samples <- nrow(data_dims)
  
  anchors <- c()
  NN_ids <- list()
  # vector of how many housing the i'th element belongs to
  housing <- rep(0,n_samples)
  
  while(any(housing < coverage)) {
    homeless <- (1:n_samples)[housing < coverage]
    if (length(homeless) == 1) {anchor <- homeless
    } else {anchor <- sample(c(homeless), 1)}
    
    anchors <- c(anchors, anchor)
    if (verbose) {
      print(paste("Anchor #",length(anchors),"Homeless:",length(homeless)))
    }
    
    if (mode=="fr") {
      # Fixed-radius mode
      NNfr <- dbscan::frNN(data_dims, eps=R, sort=FALSE, query=data_dims[anchor,])
      NN_anchor <- NNfr$id[[1]]
      #NN_ids <- NNfr$id
      # housing[NN_anchor] <- housing[NN_anchor] + 1
      # n_neighbors <- c(n_neighbors, length(NN_anchor))
    } else if (mode=="k") {
      # k-nearest mode
      NNk <- dbscan::kNN(data_dims, k=k, sort=FALSE, query=data_dims[anchor,])
      NN_anchor <- NNk$id[1,] %>% unname()
      # housing[NN_anchor] <- housing[NN_anchor] + 1
      # n_neighbors <- c(n_neighbors, k)
    } else if (mode=="hybrid") {
      # hybrid fr+k mode
      # computes NN with fixed-radius
      NNfr <- dbscan::frNN(data_dims, eps=R, sort=FALSE, query=data_dims[anchor,])
      NN_anchor <- NNfr$id[[1]]
      #NN_ids <- NNfr$id
      # computes NN with k-nearest for points with <k fixed-radius neighbors
      if (length(NN_anchor) < k) {
        NNk <- dbscan::kNN(data_dims, k=k, sort=FALSE, query=data_dims[anchor,])
        NN_anchor <- split(NNk$id,seq_len(nrow(NNk$id)))[[1]]
      }
      #n_neighbors <- c(n_neighbors, length(NN_anchor))
    }
    housing[ NN_anchor ] <- housing[ NN_anchor ] + 1
    NN_ids[[anchor]] <- NN_anchor
  }
  
  # fills out NN_ids with NULLs until n_samples is reached
  # if (length(NN_ids) != n_samples) {
  #   for (i in (length(NN_ids)+1):n_samples) { NN_ids[i] <- list(NULL) }
  # }
  # gets number of neighbors
  data$n_neighbors <- as.numeric(NA)
  data$n_neighbors[anchors] <- sapply(NN_ids, length)[anchors]
  
  print("Computing correlation within each neighborhood")
  r_values <- c()
  # computes correlation between phenotype and PGS for phenotype
  for (i in anchors) {
    r <- cor(data[ NN_ids[[i]], ][[ col_pheno ]],
             data[ NN_ids[[i]], ][[ col_PGS ]])
    r_values[i] <- r
  }
  data$locPGSacc <- as.numeric(NA)
  data$locPGSacc[anchors] <- r_values[anchors]
  
  return(data)
}
