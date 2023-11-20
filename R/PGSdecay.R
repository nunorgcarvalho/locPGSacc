#' @title PGSdecay
#' @description Returns statistics related to binned PGS accuracy as a function of some continuous/dimensional variable.
#' @inherit locPGSacc author
#' 
#' @details Takes a dataset with columns for a phenotype and its PGS as well as a 
#' column for a dimensional variable of interest (e.g. genetic distance) and splits the
#' sample into equal-length bins. PGS accuracy is calculated within each bin
#' and the rate of change of accuracy across the dimensional variable is computed.
#' 
#' @inheritParams get_NNs
#' @inheritParams get_accuracy
#' @inheritParams bin_dim
#' @param col_dim character: column name of the dimensional variable for which to compute bins from
#' @param return_objects (optional) logical: whether the function should return the outputs of [cor.test()] and [lm()] directly, rather than extracting the most important metrics 
#' 
#' @return Returns a nested list with statistics related to PGS decay across bins:
#' \itemize{
#'   \item bin_data: tibble containing data for each distance bin:
#'    \itemize{
#'      \item bin = range of distance bin
#'      \item min = minimum distance within bin
#'      \item max = max distance within bin
#'      \item N = number of samples in bin
#'      \item r = correlation between phenotype and PGS within bin
#'      \item r_lower = lower r estimate from 95% confidence interval
#'      \item r_upper = upper r estimate from 95% confidence interval
#'      \item median = median value of dimensional variable of bin's samples
#'    }
#'   \item lm: computes a linear regression for bin_R2 ~ dim_dist across all bins
#'    \itemize{
#'      \item intercept = intercept of line of best fit
#'      \item m = slope of line of best fit
#'      \item m_se = standard error of the slope of line of best fit
#'      \item p = p-value for slope
#'      \item m_hat = standardized slope; m multiplied by the windowed range of 
#'            the dimensional variable and divided by the cor(pheno, PGS) among
#'            reference population (bin with highest r_lower)
#'      \item m_hat_se = standard error of the standardized slope m_hat
#'    }
#'   \item global: computes cor(pheno, PGS) on all samples
#'    \itemize{
#'      \item r = Pearson's correlation coefficient
#'      \item p = p-value
#'    }
#'   
#' } 
#' 
#'
#' @export
#' 
#' @import tidyverse
#' @import psychometric

PGSdecay <- function (
    data,
    col_dim,
    col_pheno,
    col_PGS,
    i_omit = c(),
    bins = 15,
    ref_window = 0.95,
    return_objects = FALSE
) {
  
  # output list is established ####
  output <- list(
    bin_data = tibble(),
    lm = list(),
    global = list()
  )
  
  # gets list of samples with missing data
  i_omit <- (1:nrow(data))[is.na(data[[col_pheno]]) | is.na(data[[col_PGS]]) | is.na(data[[col_dim]])]
  if (length(i_omit) > 0) {data <- data[-i_omit,]}
  # removes missing data and renames columns of interest
  data <- data %>% dplyr::select(dim = !!sym(col_dim),
                                    pheno = !!sym(col_pheno),
                                    PGS = !!sym(col_PGS))
  # splits individuals into groups based on dim variable
  data$dim_group <- bin_dim(data, bins = bins, ref_window = ref_window)
  
  # gets per-bin performance
  bin_data <- get_bin_data(data)
  
  # computes linear regression of bin PGS accuracy against the median bin distance
  # uses sample size of bin as weight
  lm1 <- lm(R2 ~ median, data = bin_data, weights = bin_data$N)
  m <- lm1$coefficients[[2]]
  m_se <- summary(lm1)$coefficients[2,2]
  
  # gets standardized m and m_se
  m_hat.list <- standardize_m(bin_data, m, m_se)
  
  # computes correlation between phenotype and PGS for entire sample
  cor1 <- cor.test(data$pheno, data$PGS)
  
  # saves data to output list
  output$bin_data <- bin_data
  if (!return_objects) {
    output$lm$intercept<- lm1$coefficients[[1]]
    output$lm$m <- m
    output$lm$m_se <- m_se
    output$lm$p <- summary(lm1)$coefficients[2,4]
  } else {output$lm <- lm1}
  output$lm$m_hat <- m_hat.list$m_hat
  output$lm$m_hat_se <- m_hat.list$m_hat_se
  if (!return_objects) {
    output$global$r <- cor1$estimate
    output$global$R2 <- cor1$estimate^2
    output$global$p <- cor1$p.value
  } else {output$global <- cor1}
  
  
  return(output)
  
}

#' @title bin_dim
#' @description splits a complete, renamed dataset by bins along some variable
#' @inherit locPGSacc author
#' 
#' @details Takes a dataset with complete cases and with column names standardized
#' ('dim','pheno','PGS') and splits individuals into bins depending on their 'dim'
#' variable. Internal function.
#' 
#' @inheritParams PGSdecay
#' @param bins (optional) integer: number of bins to split up the middle of the sample by according to the dimensional variable. The leftmost and rightmost bins will each contain N/bins individuals. The remaining bins will be of varying sample size but equal range  
#' @param ref_window (optional) numeric: proportion of the distribution, centered at the median, to include consider when making the non-leftmost and -rightmost bins. Also affects standardization of portability slope
#' 
#' @return Returns a vector of dim groups to apply directly to the dataset:
#' 
#' @import tidyverse

bin_dim <- function(data,
                    bins = 15,
                    ref_window = 0.95
) {
  
  # splits up dimension variable according to bins
  range. <- quantile(data$dim, c(0.5 - ref_window/2, 0.5 + ref_window/2))
  breaks <- c( min(data$dim), seq(range.[1], range.[2], diff(range.)/(bins-2)), max(data$dim) )
  # assigns individuals into bins
  dim_group <- cut(data$dim, breaks = breaks, include.lowest = TRUE)
  
  return(dim_group)
}

#' @title get_bin_data
#' @description makes the bin_data tibble used in other functions
#' @inherit locPGSacc author
#' 
#' @details Takes a dataset with complete cases and with column names standardized
#' ('dim','pheno','PGS', 'dim_group') and computes per-bin statistics
#' 
#' @inheritParams PGSdecay
#'  
#' @return Returns a tibble with an observation for each bin
#' 
#' @import tidyverse

get_bin_data <- function(data
) {
  # creates empty tibble for later data
  bin_data <- tibble(bin = as.character(),
                     min = as.numeric(),
                     max = as.numeric(),
                     N = as.numeric(),
                     r = as.numeric(),
                     r_lower = as.numeric(),
                     r_upper = as.numeric(),
                     R2 = as.numeric(),
                     R2_lower = as.numeric(),
                     R2_upper = as.numeric() )
  # loops through bins ####
  for (bin in levels(data$dim_group)) {
    data_bin <- data %>% filter(dim_group == bin)
    # skips correlation in bin if less than 2 samples are present in bin
    if (nrow(data_bin) < 3) {
      cor1 <- list("estimate" = as.numeric(NA),
                   "conf.int" = as.numeric(c(NA,NA)) )
      R2 <- as.numeric(NA)
      R2_CI <- list("LCL"=as.numeric(NA),"UCL"=as.numeric(NA))
    } else {
      # computes correlation between phenotype and PGS for the bin
      cor1 <- cor.test(data_bin$PGS,data_bin$pheno)
      # computes approximate R2 confidence interval
      R2 <- cor1$estimate^2
      R2_CI <- CI.Rsq(R2, n = cor1$parameter+1, k=1, level=0.95)
    }
    # adds data to bin_data tibble
    bin_data <- bin_data %>% add_row(
      bin = bin,
      min = min(data_bin$dim),
      max = max(data_bin$dim),
      N = nrow(data_bin),
      r = cor1$estimate,
      r_lower = cor1$conf.int[1],
      r_upper = cor1$conf.int[2],
      R2 = R2,
      R2_lower = R2_CI[["LCL"]],
      R2_upper = R2_CI[["UCL"]]
    )
  }
  # computes median dimension variable value within each bin
  bin_data$median <- as.numeric(NA)
  bin_data$median[!is.na(bin_data$r)] <- ( data %>% group_by(dim_group) %>% summarize(median = median(dim)) )$median
  
  return(bin_data)
}

#' @title standardize_m
#' @description standardizes the portability slope
#' @inherit locPGSacc author
#' 
#' @details standardized slope; m multiplied by the windowed range of 
#'            the dimensional variable and divided by the cor(pheno, PGS) among
#'            reference population (bin closest to distance 0). This function also
#'            works for the standard error of the slope if provided
#' 
#' @param bin_data tibble: the bin_data tibble produced by [get_bin_data()]
#' @param m numeric: the unstandardized portability slope
#' @param m_se (optional) numeric: the unstandardized portability slope's standard error
#' 
#'  
#' @return Returns a list with the standardized portability slope (m_hat) and
#' the standardized portability slope's standard error (m_hat_se) if provided
#' 
#' @import tidyverse

standardize_m <- function(bin_data,
                          m,
                          m_se = NA
) {
  bins <- nrow(bin_data)
  # gets PGS accuracy of bin closest to zero (by median)
  i_ref <- order(bin_data$r_lower, decreasing=TRUE)
  #R2_ref <- bin_data$R2[i_ref[1]]
  R2_ref <- weighted.mean(bin_data$R2[i_ref[1:ceiling(bins/5)]],
                          bin_data$N[i_ref[1:ceiling(bins/5)]])
  # adjust m slope by the r_ref and range
  range.med <- c(bin_data$min[2],bin_data$max[bins-1])
  m_hat <- m * diff(range.med) / R2_ref %>% unname()
  # saves to output list
  out <- list(m_hat = m_hat)
  
  if (!is.na(m_se)) {
    m_hat_se <- m_se * diff(range.) / R2_ref %>% unname()
    out$m_hat_se <- m_hat_se
  }
  
  return(out)
} 