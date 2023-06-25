# small helper function that takes an NN_ids object and adds element i to the
# i'th list in the object. Returns NN_ids object.

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