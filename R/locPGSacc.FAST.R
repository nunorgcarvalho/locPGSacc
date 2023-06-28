#' @title locPGSacc.FAST
#' @description Fast version of [locPGSacc()] function for large N datasets
#' @inherit locPGSacc author
#' 
#' @details For large sample sizes (N>20,000), the full [locPGSacc()] function
#' may take too long to run and require too much memory. This function solves this
#' issue by only building neighborhoods for much smaller amount of points called
#' anchors. This function runs an algorithm to select the optimal set of anchors
#' to build neighborhoods for such that it is representative of the entire
#' dimensional space of the data set. Right now only works with quantitative
#' traits.
#' 
#' @inheritParams locPGSacc
#' @param multiplier integer: number of times each point needs to be assigned to a neighborhood before algorithm ends. Bigger = slower = more anchors.
#' @param coverage numeric: proportion of points that need to meet the multiplier parameter before algorithm ends. Bigger = slower = more anchors
#' @param seed numeric: seed used for algorithm
#' @param verbose logical: whether to print out more detailed algorithm-related information
#' 
#' @inherit locPGSacc return
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
  # separates samples with missing phenotype or PGS data
  i_NA <- which(is.na(data[[col_pheno]]) | is.na(data[[col_PGS]]))
  data_NA <- data[i_NA,] %>% select(-any_of(cols_needed))
  data <- data[-i_NA,]
  if (nrow(data_NA) > 0 ) {warning(paste0(nrow(data_NA), " sample(s) have missing PGS or phenotype data and won't be considered in algorithm."))}
  
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
    # once again checks for missing NA values, but should already be covered earlier
    data_cor <- data[ NN_ids[[i]], c(col_pheno, col_PGS)] %>% drop_na()
    r <- cor(data_cor[[ col_pheno ]],
             data_cor[[ col_PGS ]])
    r_values[i] <- r
  }
  data$locPGSacc <- as.numeric(NA)
  data$locPGSacc[anchors] <- r_values[anchors]
  
  # reappends missing data rows at the end
  if (nrow(data_NA)>0) {data <- data %>% add_row(data_NA, n_neighbors=NA, locPGSacc=NA)}
  
  return(data)
}
