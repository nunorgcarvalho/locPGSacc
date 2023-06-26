#' @title add_self2neighborhood
#' @description Takes an NN_ids object from [locPGSacc()] and adds element i to the i'th NN vector
#' @inherit locPGSacc author
#' 
#' @param NN_ids NN_ids object from locPGSacc()
#' @param i_skip (optional) integer vector: list of indices to skip when adding element i to its own NN vector
#' 
#' @returns Returns inputted NN_ids object but with element i added to the i'th NN vector

add_self2neighborhood <- function(NN_ids, i_skip=c()) {
  # skips denoted list
  i2check <- (1:length(NN_ids))[!(1:length(NN_ids) %in% i_skip)]
  for (i in i2check) {
    # only adds i to neighborhood if not already present for some reason
    if (!any(NN_ids[[i]] == i)) {
      NN_ids[[i]] <- c(i, NN_ids[[i]])
    }
  }
  return(NN_ids)
}