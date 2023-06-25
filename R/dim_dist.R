# returns inputted data table with 'dim_dist' column

library(tidyverse)

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