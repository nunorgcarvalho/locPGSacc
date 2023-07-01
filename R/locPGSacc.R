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
#' of neighbors of each point appended as columns. For smaller sample sizes
#' (N < 20,000), the method='FULL' algorithm is likely fine, but larger sample
#' sizes become too slow and resource-intensive to deal with; use method='FAST'
#' instead.
#' 
#' @inheritParams get_NNs
#' @inheritParams get_accuracy
#' @param col_PGSacc character: column name of the outputted local PGS accuracy
#' @param NN_ids (optional) nested list outputted from [get_NNs()]; allows for reuse of neighborhoods
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
    col_pheno, # for r
    col_PGS, # for r
    col_PGSacc = "locPGSacc", # for r
    method = "FAST",
    R = -1,
    k = -1,
    mode = "hybrid",
    i_omit = c(),
    NN_ids = list(),
    force_FULL = FALSE,
    # get_NNs.FAST specific parameters
    multiplier = 1,
    coverage = 1,
    seed = NA,
    verbose = FALSE
) {
  # Checks if data table already contains an output column
  cols_needed <- c("n_neighbors",col_PGSacc)
  if ( any(cols_needed %in% colnames(data)) ) {
    cols_duplicated <- cols_needed[cols_needed %in% colnames(data)]
    warning(paste0(c("Data table already contains at least one column used for output purposes:", cols_duplicated), collapse=" "))
  }
  
  if (length(NN_ids) == 0) {
    # gets indices of samples with missing phenotype of PGS data
    i_omit <- (1:nrow(data))[is.na(data[[col_pheno]]) | is.na(data[[col_PGS]])]
    
    # uses get_NNs to retrieve nearest neighbors
    NN_ids <- get_NNs(
      data,
      col_dims = col_dims,
      method = method,
      R = R,
      k = k,
      mode = mode,
      i_omit = i_omit,
      force_FULL = force_FULL,
      multiplier = multiplier,
      coverage = coverage,
      seed = seed,
      verbose = verbose
    )
  } else {
    # throws error if inputted NN_ids does not match data
    if (length(NN_ids) != nrow(data)) {stop("Length of NN_ids does not match number of rows in data")}
  }
  data$n_neighbors <- sapply(NN_ids, length)
  data$n_neighbors[data$n_neighbors==0] <- NA
  
  # computes correlation between phenotype and PGS for phenotype
  print("Computing correlation within each neighborhood")
  r_values <- c()
  for (i in 1:length(NN_ids)) {
    NN_ids_i <- NN_ids[[i]]
    # skips if no NNs were computed
    if (is.null(NN_ids_i)) {next}
    # once again checks for missing NA values, but should already be covered earlier
    # data_cor <- data[ NN_ids_i, c(col_pheno, col_PGS)] %>% drop_na()
    # r <- cor(data_cor[[ col_pheno ]],
    #          data_cor[[ col_PGS ]])
    r <- get_accuracy(data, NN_ids_i, col_pheno=col_pheno, col_PGS=col_PGS)
    r_values[i] <- r
  }
  data[[col_PGSacc]] <- r_values
  
  # returns full data set
  return(data)
}




# accuracy function ####


#' @title get_accuracy
#' @description Computes correlation between phenotype and PGS for given neighborhood
#' @inherit locPGSacc author
#' 
#' @details This function takes a dataset containing columns for a phenotype and
#' a PGS for the phenotype as well as a list of indices in the dataset (composing
#' the "neighborhood") for which to compute the correlation between phenotype
#' and PGS for.
#' 
#' @inheritParams locPGSacc
#' @param NN_ids_i: numeric vector: indices in 'data' that are part of neighborhood to compute correlation within
#' @param col_pheno character: column name of the phenotype of interest
#' @param col_PGS character: column name of the polygenic scores for the phenotype of interest
#' 
#' @return Returns correlation value between phenotype and PGS within neighborhood
#' 
#' @import tidyverse
get_accuracy <- function(
    data,
    NN_ids_i,
    col_pheno,
    col_PGS
) {
  # removes potential NA values (this should have already been taken care of)
  data_cor <- data[ NN_ids_i, c(col_pheno, col_PGS)] %>% drop_na()
  r <- cor(data_cor[[ col_pheno ]],
           data_cor[[ col_PGS ]])
  
  # returns r
  return(r)
}