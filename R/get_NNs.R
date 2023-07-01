#' @title get_NNs
#' @description Returns nested list of samples' nearest neighbors
#' @author Nuno R. G. Carvalho: \email{nunocarvalho@@gatech.edu}
#' 
#' @details This function takes a dataset containing columns for a phenotype,
#' a PGS for the phenotype, and coordinates for some dimensional-space (e.g.
#' genetic principal component space) and computes neighborhoods. A point's
#' neighborhood can be defined in different ways, but broadly consists of points
#' close to the original point chosen (called the anchor point).
#' 
#' @param data data table containing all necessary columns and rows
#' @param col_dims character vector: name of the dimension column(s)
#' @param R (optional for mode='k') numeric: radius from each point to build a neighborhood from. Only needed if using 'fr' or 'hybrid' mode
#' @param k (optional for mode='fr') integer: number of closest neighbors (including self) to build neighborhood from. Only needed if using 'k' or 'hybrid' mode
#' @param mode (optional) character: mode of building neighborhoods. One of:
#' \itemize{
#'   \item 'fr' = fixed-radius mode. Each anchor's neighborhood is composed of every point within a radius R of the anchor point.
#'   \item 'k' = k-nearest mode. Each anchor's neighborhood is composed of the closest k points to the anchor point.
#'   \item 'hybrid' = fr+k mode. Runs fixed-radius mode and then runs k-nearest mode on neighborhoods with sizes less than k. Recommended for reducing noise.
#' }
#' @param i_omit (optional) numeric vector: indices of data rows for which to not consider at all in nearest neighbor algorithm
#' @param force_large_N (optional) logical: whether to continue with function if data's sample size is large (N>20,000). It is recommended to use [locPGSacc.FAST()] in such situations.
#' 
#' @return Returns a nested list, where the i'th element contains a numeric vector of
#' the i'th point's neighbors (including element i itself)
#' 
#'
#' @export
#' 
#' @import tidyverse
#' @import dbscan

get_NNs <- function(
    data,
    col_dims,
    R = NA,
    k = NA,
    mode = "hybrid",
    i_omit = c(),
    force_large_N = FALSE
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
  # stops if N is large and user didn't force slow algorithm
  if (n_samples > 20000 & !force_large_N) {
    stop("N is (probably) too big for this algorithm. Use locPGSacc.FAST() instead,
         or set 'force_largeN'to TRUE if you are sure you want to proceed")
  }
  # warns if k > n_samples
  if (k >= n_samples) {warning("k=",k," is bigger than number of non-omitted samples")}
  
  print(paste0("Computing nearest neighbors for ",n_samples," samples using ",mode,
               "-mode with the settings",
               ifelse(mode!="k",paste(" R =",R),""),
               ifelse(mode!="fr",paste(" k =",k),"")))
  # actually computes NN_ids
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
  
  # maps post-omit indices back to pre-omit indices
  NN_ids_keep <- sapply(NN_ids_keep, function(x) i_keep[x])
  
  # makes NN_ids for all samples in full data, setting omitted samples to NULL
  NN_ids <- vector("list", nrow(data))
  NN_ids[i_keep] <- NN_ids_keep
  
  # returns NN_ids2
  return(NN_ids)
}