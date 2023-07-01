#' @title get_NNs
#' @description Returns nested list of samples' nearest neighbors
#' @inherit locPGSacc author
#' 
#' @details This function takes a dataset containing columns for a phenotype,
#' a PGS for the phenotype, and coordinates for some dimensional-space (e.g.
#' genetic principal component space) and computes neighborhoods. A point's
#' neighborhood can be defined in different ways, but broadly consists of points
#' close to the original point chosen (called the anchor point).
#' 
#' @inheritParams get_NNs.FULL
#' @inheritParams get_NNs.FAST
#' @param data data table containing all necessary columns and rows
#' @param col_dims character vector: name of the dimension column(s)
#' @param method character vector: method of sampling points to build neighborhoods for. One of:
#' \itemize{
#'   \item 'FULL' = builds neighborhoods for every points in the dataset (with non-NA data)
#'   \item 'FAST' = builds neighborhoods for a much smaller but representative sample of the data's dimensional space. Algorithm can be tweaked further (see below)
#' }
#' @param R (optional for mode='k') numeric: radius from each point to build a neighborhood from. Only needed if using 'fr' or 'hybrid' mode
#' @param k (optional for mode='fr') integer: number of closest neighbors (including self) to build neighborhood from. Only needed if using 'k' or 'hybrid' mode
#' @param mode (optional) character: mode of building neighborhoods. One of:
#' \itemize{
#'   \item 'fr' = fixed-radius mode. Each anchor's neighborhood is composed of every point within a radius R of the anchor point.
#'   \item 'k' = k-nearest mode. Each anchor's neighborhood is composed of the closest k points to the anchor point.
#'   \item 'hybrid' = fr+k mode. Runs fixed-radius mode and then runs k-nearest mode on neighborhoods with sizes less than k. Recommended for reducing noise.
#' }
#' @param i_omit (optional) numeric vector: indices of data rows for which to not consider at all in nearest neighbor algorithm
#' @param force_FULL (optional, method='FULL') logical: whether to continue with function if data's sample size is large (N>20,000). It is recommended to use [locPGSacc.FAST()] in such situations.
#' 
#' @return Returns a nested list, where the i'th element contains a numeric vector of
#' the i'th point's neighbors (including element i itself) if computed
#' 
#'
#' @export
#' 
#' @import tidyverse
#' @import dbscan

get_NNs <- function(
    # general parameters
    data,
    col_dims,
    method = "FULL",
    R = NA,
    k = NA,
    mode = "hybrid",
    i_omit = c(),
    force_FULL = FALSE,
    # get_NNs.FAST specific parameters
    multiplier = 1,
    coverage = 1,
    seed = NA,
    verbose = FALSE
) {
  
  # Checks if correct parameter was supplied
  if (mode != "k" & is.na(R)) {stop("R must be supplied")}
  if (mode != "fr" & is.na(k)) {stop("k must be supplied")}
  
  # Checks if col_dims are in data table
  if ( !all(col_dims %in% colnames(data))) {stop("Data table does not contain all specified col_dims")}
  
  # adds samples with NA data to list of indices to omit
  i_NA <- (1:nrow(data))[rowSums(is.na(data[,col_dims]))>0]
  if (length(i_NA) > 0 ) {warning(paste0(length(i_NA), " sample(s) have missing dimensional data and won't be considered in algorithm."))}
  i_omit <- c(i_omit,i_NA) %>% unique() %>% sort()
  if (length(i_omit) == 0) {i_keep <- 1:nrow(data)
  } else {i_keep <- (1:nrow(data))[-i_omit]}
  
  # makes data table with just col_dims and no missing data
  data_dims <- data[i_keep,] %>% select(all_of(col_dims))
  n_samples <- nrow(data_dims)
  
  # warns if k > n_samples
  if (k >= n_samples) {warning("k=",k," is bigger than number of non-omitted samples")}
  
  print(paste0("Computing '",method,"' nearest neighbors for ",n_samples," samples using ",
               mode, "-mode with the settings",
               ifelse(mode!="k",paste(" R =",R),""),
               ifelse(mode!="fr",paste(" k =",k),"")))
  
  # actually computes NN_ids
  if (method %in% c("FULL","full")) {
    # stops if N is large and user didn't force slow algorithm
    if (n_samples > 20000 & !force_FULL) {
      stop("N is (probably) too big for this algorithm. Use locPGSacc.FAST() instead,
         or set 'force_largeN'to TRUE if you are sure you want to proceed")
    }
    NN_ids_keep <- get_NNs.FULL(data_dims, R=R, k=k, mode=mode)
  } else if (method %in% c("FAST","fast")) {
    NN_ids_keep <- get_NNs.FAST(data_dims, R=R, k=k, mode=mode,
                                multiplier=multiplier, coverage=coverage,
                                seed=seed, verbose=verbose)
  } else { stop("Missing correct 'method' of NN computation") }
  
  # maps post-omit indices back to pre-omit indices
  NN_ids_keep <- sapply(NN_ids_keep, function(x) i_keep[x])
  
  # makes NN_ids for all samples in full data, setting omitted samples to NULL
  NN_ids <- vector("list", nrow(data))
  NN_ids[i_keep] <- NN_ids_keep
  
  # returns NN_ids
  return(NN_ids)
}




# FULL NN function ####

#' @title get_NNs.FULL
#' @description NN algorithm that builds all possible neighborhoods
#' @inherit locPGSacc author
#' 
#' @details This NN algorithm builds a neighborhood for every point in the data
#' set (with non-NA data)
#' 
#' @inheritParams get_NNs
#' @param data_dims data table containing all necessary columns and rows. MUST HAVE NO NA VALUES
#'  
#' @inherit get_NNs return
#' 
#' @import tidyverse
#' @import dbscan

get_NNs.FULL <- function(
    data_dims, # can't contain any missing data
    R = NA,
    k = NA,
    mode = "hybrid"
) {
  if (mode=="fr") {
    # Fixed-radius mode
    NNfr <- dbscan::frNN(data_dims, eps=R, sort=FALSE)
    NN_ids_keep <- NNfr$id %>% add_self2neighborhood()
  } else if (mode=="k") {
    # k-nearest mode
    NNk <- dbscan::kNN(data_dims, k=k-1, sort=FALSE)
    NN_ids_keep <- split(NNk$id,seq_len(nrow(NNk$id))) %>% add_self2neighborhood()
  } else if (mode=="hybrid") {
    # hybrid fr+k mode
    # computes NN with fixed-radius
    NNfr <- dbscan::frNN(data_dims, eps=R, sort=FALSE)
    NN_ids_keep <- NNfr$id
    # computes NN with k-nearest for points with <k fixed-radius neighbors
    small_neighborhoods <- (1:nrow(data_dims))[sapply(NNfr$id, length) < k]
    if (length(small_neighborhoods > 0)) {
      NNk <- dbscan::kNN(data_dims, k=k, sort=FALSE, query=data_dims[small_neighborhoods,])
      NN_ids_keep[small_neighborhoods] <- split(NNk$id,seq_len(nrow(NNk$id)))
    }
    NN_ids_keep <- add_self2neighborhood(NN_ids_keep)
  }
  
  # returns NN_ids_keep
  return(NN_ids_keep)
}


#' @title get_NNs.FAST
#' @description NN algorithm that builds all possible neighborhoods
#' @inherit locPGSacc author
#' 
#' @details This NN algorithm builds a neighborhood for every point in the data
#' set (with non-NA data)
#' 
#' @inheritParams get_NNs
#' @param data_dims data table containing all necessary columns and rows. MUST HAVE NO NA VALUES
#' @param multiplier (optional, method='FAST') integer: number of times each point needs to be assigned to a neighborhood before algorithm ends. Bigger = slower = more anchors.
#' @param coverage (optional, method='FAST') numeric: proportion of points that need to meet the multiplier parameter before algorithm ends. Bigger = slower = more anchors
#' @param seed (optional, method='FAST') numeric: seed used for algorithm
#' @param verbose (optional, method='FAST') logical: whether to print out more detailed algorithm-related information
#'  
#' @inherit get_NNs return
#' 
#' @import tidyverse
#' @import dbscan

# FAST NN function ####
get_NNs.FAST <- function(
    data_dims, # can't contain any missing data
    R = NA,
    k = NA,
    mode = "hybrid",
    multiplier = 1,
    coverage = 1,
    seed = NA,
    verbose = FALSE
) {
  n_samples <- nrow(data_dims)
  # empty vectors for later loop
  anchors <- c()
  NN_ids_keep <- vector("list", n_samples)
  housing <- rep(0,n_samples) # vector of how many housing the i'th element belongs to
  
  if (verbose) {print(paste0("Looping until less than ",ceiling(n_samples*(1-coverage))," unhoused remain"))}
  # actual loop
  while(sum(housing < multiplier) > n_samples*(1-coverage)) {
    unhoused <- (1:n_samples)[housing < multiplier]
    if (length(unhoused) == 1) {anchor <- unhoused
    } else {anchor <- sample(c(unhoused), 1)}
    
    anchors <- c(anchors, anchor)
    if (verbose) { print(paste("Anchor #",length(anchors),"Unhoused:",length(unhoused))) }
    
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
    NN_ids_keep[[anchor]] <- NN_anchor
  }
  print(paste0("Formed ", length(anchors), " neighborhoods with a coverage of ",
               coverage," and multiplier of ", multiplier))
  
  return(NN_ids_keep)
}