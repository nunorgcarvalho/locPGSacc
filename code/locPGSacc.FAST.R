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
    mode = "fr", # fr=fixed-radius, k=k closest neighbors, hybrid=highest n_neighbors of fr and k
    coverage = 1,
    seed = NA
) {
  if (!is.na(seed)) { set.seed(seed) }
  
  data_dims <- data %>% select(any_of(col_dims))
  n_samples <- nrow(data_dims)
  
  anchors <- c()
  n_neighbors <- c()
  # vector of how many housing the i'th element belongs to
  housing <- rep(0,n_samples)
  
  while(any(housing < coverage)) {
    homeless <- (1:n_samples)[housing < coverage]
    if (length(homeless) == 1) {anchor <- homeless
    } else {anchor <- sample(c(homeless), 1)}
    
    anchors <- c(anchors, anchor)
    print(paste("Anchor #",length(anchors),"Homeless:",length(homeless)))
    
    # NNfr <- dbscan::frNN(data_dims, eps=R, sort=FALSE, query=data_dims[anchor,])
    # housing[NNfr$id[[1]]] <- housing[NNfr$id[[1]]] + 1
    # n_neighbors <- c(n_neighbors, length(NNfr$id[[1]]))
    
    if (mode=="fr") {
      # Fixed-radius mode
      #print(paste0("Computing fixed-radius=",R," nearest neighbors for ",n_samples," samples"))
      NNfr <- dbscan::frNN(data_dims, eps=R, sort=FALSE, query=data_dims[anchor,])
      NN_ids <- NNfr$id
      housing[NN_ids[[1]]] <- housing[NN_ids[[1]]] + 1
      n_neighbors <- c(n_neighbors, length(NN_ids[[1]]))
    } else if (mode=="k") {
      # k-nearest mode
      #print(paste0("Computing k=",k," nearest neighbors for ",n_samples," samples"))
      NNk <- dbscan::kNN(data_dims, k=k, sort=FALSE, query=data_dims[anchor,])
      NN_ids <- split(NNk$id,seq_len(nrow(NNk$id)))
      housing[NNfr$id[[1]]] <- housing[NNfr$id[[1]]] + 1
      n_neighbors <- c(n_neighbors, k)
    } else if (mode=="hybrid") {
      # hybrid fr+k mode
      #print(paste0("Computing hybrid fixed-radius=",R," & k=",k," nearest neighbors for ",n_samples," samples"))
      # computes NN with fixed-radius
      NNfr <- dbscan::frNN(data_dims, eps=R, sort=FALSE, query=data_dims[anchor,])
      NN_ids <- NNfr$id
      # computes NN with k-nearest for points with <k fixed-radius neighbors
      if (length(NN_ids[[1]]) < k) {
        NNk <- dbscan::kNN(data_dims, k=k, sort=FALSE, query=data_dims[anchor,])
        NN_ids[1] <- split(NNk$id,seq_len(nrow(NNk$id)))
        #housing[NNk$id[[1]]] <- housing[NNk$id[[1]]] + 1
      }
      housing[ NN_ids[[1]] ] <- housing[ NN_ids[[1]] ] + 1
      n_neighbors <- c(n_neighbors, length(NN_ids[[1]]))
    }
    
  }
}