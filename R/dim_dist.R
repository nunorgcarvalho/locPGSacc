#' @title dim_dist
#' @description Computes each sample's distance from a specified reference point (centroid if nto given)
#' @inherit locPGSacc author
#' 
#' @details This function takes a dataset containing sample coordinates in some
#' dimensional space and computes each sample's distance from a reference point.
#' If not specified, the reference point will be automatically set to the
#' centroid (mean position) of all the samples. If specified, the centroid will
#' only be computed for a subset of samples.
#' 
#' @param data data table containing all necessary columns and rows
#' @param col_dims character vector: name of the dimension column(s)
#' @param reference_point (optional) specifies what reference point is used to calculate sample distance from. One of:
#' \itemize{
#'   \item default: uses centroid (mean position) of all samples. Specify which samples to consider for centroid using 'centroid_indices' parameter.
#'   \item 0: origin of the space
#'   \item numeric vector (same length as col_dims)
#' }
#' @param centroid_indices (optional) integer or logical vector: indices/mask of the data's samples to use for centroid calculation. All samples are used by default
#' @param col_dist (optional) character: output column name of sample's distance in space
#' 
#' @export
#' 
#' @import tidyverse

dim_dist <- function(
    data,
    col_dims,
    reference_point = NA, # must be same length as number of dimensions, except 0 = origin
    centroid_indices = NA, # vector of indices to compute centroid for (if no reference point given). Leave NA for all points
    col_dist = "dim_dist"
) {
  # Checks if col_dims are in data table
  if ( !all(col_dims %in% colnames(data))) {stop("Data table does not contain all specified col_dims")}
  # Checks if data table already contains an output column
  cols_needed <- c("n_neighbors","locPGSacc")
  if ( col_dist %in% colnames(data) ) {
    warning(paste0("Data table already contains column used for output: ", col_dist))
  }
  
  data_dims <- data %>% select(any_of(col_dims))
  
  # if not supplied, reference point is set to mean of data (centroid)
  if (is.na(reference_point)) {
    print("Setting reference point to the centroid of specified points (mean of dimensions)")
    if (all(!is.na(centroid_indices))) {
      reference_point <- colMeans(data_dims[centroid_indices,])
    } else {
      reference_point <- colMeans(data_dims)
    }
  }
  # reports error if length of reference_point vector is not the same as number of dimensions
  if ( (length(reference_point) != length(col_dims)) & (!identical(reference_point, 0))) {
    stop(paste0("Number of dimensions (",length(col_dims),") does not match dimensions of reference point (",length(reference_point),")"))
  }
  
  # Computes distance from reference point
  print(paste0("Computing distance from specified reference point"))
  data[,col_dist] <- sqrt(rowSums((data_dims - reference_point)^2))
  
  return (data)
}