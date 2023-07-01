#' @title locPGSacc
#' @description Computes local PGS accuracy for all points in inputed data set
#' @author Nuno R. G. Carvalho: \email{nunocarvalho@@gatech.edu}
#' 
#' @details This function takes a dataset containing columns for a phenotype,
#' a PGS for the phenotype, and coordinates for some dimensional-space (e.g.
#' genetic principal component space) and computes the correlation between the 
#' phenotype and PGS within a 'neighborhood' centered at each point. This is called
#' local PGS accuracy. A point's neighborhood can be defined in different ways,
#' but broadly consists of points close to the original point chosen. This
#' function returns the inputted data set with the local PGS accuracy and number
#' of neighbors of each point appended as columns. For large sample sizes 
#' (N > 20,000), it is recommended to use [locPGSacc.FAST()] instead
#' 
#' @inheritParams get_NNs
#' @param col_pheno character: column name of the phenotype of interest
#' @param col_PGS character: column name of the polygenic scores for the phenotype of interest
#' @param col_PGSacc character: column name of the outputted local PGS accuracy
#' 
#' @return Returns inputted 'data' table but with the following columns appended:
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

locPGSacc <- function (
    data,
    col_dims,
    col_pheno,
    col_PGS,
    col_PGSacc = "locPGSacc",
    R = -1,
    k = -1,
    mode = "hybrid", # fr=fixed-radius, k=k closest neighbors, hybrid=highest n_neighbors of fr and k
    force_large_N = FALSE
) {
  # Checks if data table already contains an output column
  cols_needed <- c("n_neighbors",col_PGSacc)
  if ( any(cols_needed %in% colnames(data)) ) {
    cols_duplicated <- cols_needed[cols_needed %in% colnames(data)]
    warning(paste0(c("Data table already contains at least one column used for output purposes:", cols_duplicated), collapse=" "))
  }
  
  # uses get_NNs to retrieve nearest neighbors
  NN_ids <- get_NNs(
    data,
    col_dims=col_dims,
    R = R,
    k = k,
    mode = mode,
    i_omit = c(),
    force_large_N = force_large_N
  )
  data$n_neighbors <- sapply(NN_ids, length)
  # adds lis
  
  # # Checks if correct parameter was supplied
  # if (mode != "k" & R==-1) {stop("R must be supplied")}
  # if (mode != "fr" & k==-1) {stop("k must be supplied")}
  # 
  # # Checks if col_dims are in data table
  # if ( !all(col_dims %in% colnames(data))) {stop("Data table does not contain all specified col_dims")}
  
  
  # separates samples with missing phenotype or PGS data
  # i_NA <- which(is.na(data[[col_pheno]]) | is.na(data[[col_PGS]]))
  # data_NA <- data[i_NA,] %>% select(-any_of(cols_needed))
  # data <- data[-i_NA,]
  # if (nrow(data_NA) > 0 ) {warning(paste0(nrow(data_NA), " sample(s) have missing PGS or phenotype data and won't be considered in algorithm."))}
  
  # data_dims <- data %>% select(any_of(col_dims))
  # n_samples <- nrow(data_dims)
  # if (n_samples > 20000 & !force_large_N) {
  #   stop("N is (probably) too big for this algorithm. Use locPGSacc.FAST() instead,
  #        or set 'force_largeN'to TRUE if you are sure you want to proceed")
  #   }
  # if (k >= n_samples) {warning("k=",k," is bigger than number of samples")}
  
  # print(paste0("Computing nearest neighbors for ",n_samples," samples using ",mode,
  #              "-mode with the settings",
  #              ifelse(mode!="k",paste(" R =",R),""),
  #              ifelse(mode!="fr",paste(" k =",k),"")))
  # if (mode=="fr") {
  #   # Fixed-radius mode
  #   NNfr <- dbscan::frNN(data_dims, eps=R, sort=FALSE)
  #   NN_ids <- NNfr$id %>% add_self2neighborhood()
  #   data$n_neighbors <- sapply(NN_ids, length)
  # } else if (mode=="k") {
  #   # k-nearest mode
  #   NNk <- dbscan::kNN(data_dims, k=k-1, sort=FALSE)
  #   NN_ids <- split(NNk$id,seq_len(nrow(NNk$id))) %>% add_self2neighborhood()
  #   data$n_neighbors <- k
  # } else if (mode=="hybrid") {
  #   # hybrid fr+k mode
  #   # computes NN with fixed-radius
  #   NNfr <- dbscan::frNN(data_dims, eps=R, sort=FALSE)
  #   NN_ids <- NNfr$id
  #   # computes NN with k-nearest for points with <k fixed-radius neighbors
  #   small_neighborhoods <- (1:nrow(data))[sapply(NNfr$id, length) < k]
  #   if (length(small_neighborhoods > 0)) {
  #     NNk <- dbscan::kNN(data_dims, k=k, sort=FALSE, query=data_dims[small_neighborhoods,])
  #     NN_ids[small_neighborhoods] <- split(NNk$id,seq_len(nrow(NNk$id)))
  #   }
  #   NN_ids <- add_self2neighborhood(NN_ids)
  #   data$n_neighbors <- sapply(NN_ids, length)
  # }
  
  # computes correlation between phenotype and PGS for phenotype
  print("Computing correlation within each neighborhood")
  r_values <- c()
  for (i in 1:length(NN_ids)) {
    # skips if no NNs were computed
    if (is.null(NN_ids[[i]])) {next}
    # once again checks for missing NA values, but should already be covered earlier
    data_cor <- data[ NN_ids[[i]], c(col_pheno, col_PGS)] %>% drop_na()
    r <- cor(data_cor[[ col_pheno ]],
             data_cor[[ col_PGS ]])
    r_values[i] <- r
  }
  data[[col_PGSacc]] <- r_values
  
  # reappends missing data rows at the end
  #if (nrow(data_NA)>0) {data <- data %>% add_row(data_NA)}
  
  return(data)
}