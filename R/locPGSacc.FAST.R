#' @title locPGSacc.FAST
#' @description Fast version of locPGSacc() function for large N datasets
#' @author Nuno R. G. Carvalho: \email{nunocarvalho@@gatech.edu}
#' 
#' @details For large sample sizes (N>20,000), the full locPGSacc() function may
#' take too long to run and require too much memory. This function solves this
#' issue by only building neighborhoods for much smaller amount of points called
#' anchors. This function runs an algorithm to select the optimal set of anchors
#' to build neighborhoods for such that it is representative of the entire
#' dimensional space of the data set.
#' 
#' @param data Data table containing all necessary columns and rows
#' @param col_dims Vector of strings containing the name of the dimension columns
#' @param col_pheno String of the column name of the phenotype of interest
#' @param col_PGS String of the column name of the polygenic scores for the phenotype of interest
#' @param R Number denoting the radius from each point to build a neighborhood from. Only needed if using 'fr' or 'hybrid' mode
#' @param k Number denoting the number of closest neighbors (including self) to build neighborhood from. Only needed if using 'k' or 'hybrid' mode
#' @param mode Mode of building neighborhoods. One of:
#' \itemize{
#'   \item 'fr' = fixed-radius mode. Each anchor's neighborhood is composed of every point within a radius R of the anchor point.
#'   \item 'k' = k-nearest mode. Each anchor's neighborhood is composed of the closest k points to the anchor point.
#'   \item 'hybrid' = fr+k mode. Runs fixed-radius mode and then runs k-nearest mode on neighborhoods with sizes less than k. Recommended for reducing noise.
#' }
#' @param multiplier Integer of times each point needs to be assigned to a neighborhood before algorithm ends. Bigger = slower = more anchors.
#' @param coverage Number denoting proportion of points that need to meet the multiplier parameter before algorithm ends. Bigger = slower = more anchors
#' @param seed Seed used for algorithm
#' @param verbose Logical denoting whether to print out more detailed algorithm-related information
#' 
#' @returns Returns inputted 'data' table but with the following columns appended:
#' \itemize{
#'  \item 'n_neighbors' = number of neighbors (including self) for that point's neighborhood
#'  \item 'locPGSacc' = correlation between 'col_pheno' and 'col_PGS' columns within that point's neighborhood
#' } 
#' 
#'
#' @export
#' 
#' @import tidyverse
#' @import dbscan

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
  
  # Checks if correct parameter was supplied
  if (mode != "k" & R==-1) {stop("R must be supplied")}
  if (mode != "fr" & k==-1) {stop("k must be supplied")}
  
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
  
  print(paste0("Computing nearest neighbors for ",n_samples," samples using ",mode,
               "-mode with the settings",
               ifelse(mode!="k",paste(" R =",R),""),
               ifelse(mode!="fr",paste(" k =",k),"")))
  
  if (verbose) {print(paste0("Looping until less than ",ceiling(n_samples*(1-coverage))," unhoused remain"))}
  
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
  print(paste0("Formed ", length(anchors), " neighborhoods for a coverage of ",
               coverage," and multiplier of ", multiplier))
  
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
