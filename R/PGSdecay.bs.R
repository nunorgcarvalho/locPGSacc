#' @title PGSdecay.bs
#' @description Computes bin-based PGS decay through bootstrapping
#' @inherit locPGSacc author
#' 
#' @details Works the same as [PGSdecay] but employs stratified bootstrapping
#' to generate multiple estimated regression slopes. Instead of resampling over the 
#' entire sample, resampling is done within each bin to ensure equal representation
#' across the distance variable of interest. Best used for estimating the
#' significance of PGS decay along some variable. Can take a while to finish.
#' 
#' @inheritParams PGSdecay
#' @param B (optional) integer: number of bootstrap iterations to perform
#' @param verbose (optional) logical: whether bootstrapping iteration status should be reported to the console
#' 
#' @return returns a list containing parameters estimated from different bootstrap iterations
#' \itemize{
#'   \item m: unstandardized slopes
#'   \item m_se: unstandardized slopes' standard errors
#'   \item m_hat: standardized slopes
#'   \item m_hat_se: standardized slopes' standard errors
#' }
#' 
#'
#' @export
#' 
#' @import tidyverse
#' @import psychometric
#' 

PGSdecay.bs <- function (
    data,
    col_dim,
    col_pheno,
    col_PGS,
    B = 1000,
    i_omit = c(),
    ref_window = 0.95,
    bins = 15,
    return_objects = FALSE,
    verbose = TRUE
) {
  
  # used for estimating time left
  t1 <- Sys.time()
  
  # makes output list
  out.bs <- list(m = c(), m_se = c(), m_hat = c())
  
  # gets list of samples with missing data
  i_omit <- (1:nrow(data))[is.na(data[[col_pheno]]) | is.na(data[[col_PGS]]) | is.na(data[[col_dim]])]
  if (length(i_omit) > 0) {data <- data[-i_omit,]}
  # removes missing data and renames columns of interest
  data <- data %>% dplyr::select(dim = !!sym(col_dim),
                                 pheno = !!sym(col_pheno),
                                 PGS = !!sym(col_PGS))
  # splits individuals into groups based on dim variable
  data$dim_group <- bin_dim(data, ref_window = ref_window, bins = bins)
  
  for (b in 1:B) {
    
    # stratified bootstrap resampling
    data_bs <- data[0,]
    for (bin in levels(data$dim_group)) {
      bin_N <- sum(data$dim_group == bin)
      bs_indices <- sample(1:bin_N, bin_N, replace=TRUE)
      data_bs <- bind_rows(data_bs, data[data$dim_group == bin,][bs_indices,])
    }
    
    # gets bin data
    bin_data <- get_bin_data(data_bs)
    
    # computes linear regression of bin PGS accuracy against the median bin distance
    # uses sample size of bin as weight
    lm1 <- lm(R2 ~ median, data = bin_data, weights = bin_data$N)
    m <- lm1$coefficients[[2]]
    m_se <- summary(lm1)$coefficients[2,2]
    
    # gets standardized m and m_se
    m_hat.list <- standardize_m(data_bs, bin_data, m, m_se, ref_window)
    
    # saves to output list
    out.bs$m[b] <- m
    out.bs$m_se[b] <- m_se
    out.bs$m_hat[b] <- m_hat.list$m_hat
    out.bs$m_hat_se[b] <- m_hat.list$m_hat_se
    
    # prints to console
    if (verbose & b %% 10 == 0) {
      t_diff <- difftime(Sys.time(), t1, units='secs')
      t_left <- t_diff * (B - b) / b
      
      print(paste0('b = ',b, ' :: m_hat = ', round(m_hat.list$m_hat,4),
                   ' :: Time left ', round(t_left,1),' secs'))
    }
  }
  
  return(out.bs)
}