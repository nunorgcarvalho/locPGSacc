#' @title get_bandPGS_decay
#' @description Returns statistics related to banded PGS accuracy as a function of some dimensional variable.
#' @inherit locPGSacc author
#' 
#' @details Takes a dataset with columns for a phenotype and its PGS as well as a 
#' column for a dimensional variable of interest (e.g. genetic distance) and splits the
#' sample into equal-length bands/sections. PGS accuracy is calculated within each
#' band and the rate of change of accuracy across the dimensional variable is computed.
#' 
#' @inheritParams get_NNs
#' @inheritParams get_accuracy
#' @param col_dim character: column name of the dimensional variable for which to compute bands from
#' @param ref_window (optional) numeric: proportion of the distribution, centered at the median, to include when making bands
#' @param bands (optional) integer: number of bands (i.e. sections) to split up the middle of the sample by according to the dimensional variable. An extra band is made on each end to include dispersed samples.  
#' @param return_objects (optional) logical: whether the function should return the outputs of cor.test() and lm() directly, rather than extracting the most important metrics 
#' 
#' @return Returns a nested list with statistics related to banded PGS decay:
#' \itemize{
#'   \item band_data: tibble containing data for each distance band:
#'    \itemize{
#'      \item band = range of distance band
#'      \item min = minimum distance within band
#'      \item max = max distance within band
#'      \item N = number of samples in band
#'      \item r = correlation between phenotype and PGS within band
#'      \item r_lower = lower r estimate from 95% confidence interval
#'      \item r_upper = upper r estimate from 95% confidence interval
#'      \item median = median value of dimensional variable of band's samples
#'    }
#'   \item lm: computes a linear regression for bandPGSacc ~ dim_dist across the bands
#'    \itemize{
#'      \item intercept = intercept of line of best fit
#'      \item m = slope of line of best fit
#'      \item m_se = standard error of the slope of line of best fit
#'      \item p = p-value for slope
#'      \item m_hat = standardized slope; m multiplied by the windowed range of 
#'            the dimensional variable and divided by the cor(pheno, PGS) among
#'            reference population (band closest to distance 0)
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

get_bandPGS_decay <- function (
    data,
    col_dim,
    col_pheno,
    col_PGS,
    i_omit = c(),
    ref_window = 0.95,
    bands = 15,
    return_objects = FALSE
) {
  
  # output list is established ####
  output <- list(
    band_data = tibble(),
    lm = list(),
    global = list()
  )
  
  # gets list of samples with missing data
  i_omit <- (1:nrow(data))[is.na(data[[col_pheno]]) | is.na(data[[col_PGS]]) | is.na(data[[col_dim]])]
  # removes missing data and renames columns of interest
  data <- data[-i_omit,] %>% dplyr::select(dim = !!sym(col_dim),
                                    pheno = !!sym(col_pheno),
                                    PGS = !!sym(col_PGS))
  # gets range of dimension variable within specified central window of population distribution
  low_bound_i <- round((0.5 - ref_window/2)*length(data$dim))
  upp_bound_i <- round((0.5 + ref_window/2)*length(data$dim))
  range. = sort(data$dim)[low_bound_i:upp_bound_i] %>% range()
  # splits up dimension variable according to bands, plus <window and >window samples
  breaks <- c( min(data$dim), seq(range.[1], range.[2], diff(range.)/(bands-2)), max(data$dim) )
  # assigns individuals into bands
  data$dim_group <- cut(data$dim, breaks = breaks, include.lowest = TRUE)
  # creates empty tibble for later data
  band_data <- tibble(band = as.character(),
                      min = as.numeric(),
                      max = as.numeric(),
                      N = as.numeric(),
                      r = as.numeric(),
                      r_lower = as.numeric(),
                      r_upper = as.numeric(),
                      R2 = as.numeric(),
                      R2_lower = as.numeric(),
                      R2_upper = as.numeric() )
  # loops through band ####
  for (band in levels(data$dim_group)) {
    data_band <- data %>% filter(dim_group == band)
    # skips correlation in band if less than 2 samples are present in band
    if (nrow(data_band) < 3) {
      cor1 <- list("estimate" = as.numeric(NA),
                   "conf.int" = as.numeric(c(NA,NA)) )
    } else {
      # computes correlation between phenotype and PGS for the band
      cor1 <- cor.test(data_band$PGS,data_band$pheno)
    }
    # computes approximate R2 confidence interval
    R2 <- cor1$estimate^2
    R2_CI <- CI.Rsq(R2, n = cor1$parameter+1, k=1, level=0.95)
    # adds data to band_data tibble
    band_data <- band_data %>% add_row(
      band = band,
      min = min(data_band$dim),
      max = max(data_band$dim),
      N = nrow(data_band),
      r = cor1$estimate,
      r_lower = cor1$conf.int[1],
      r_upper = cor1$conf.int[2],
      R2 = R2,
      R2_lower = R2_CI[["LCL"]],
      R2_upper = R2_CI[["UCL"]]
    )
  }
  # computes median dimension variable value within each band
  band_data$median <- as.numeric(NA)
  band_data$median[!is.na(band_data$r)] <- ( data %>% group_by(dim_group) %>% summarize(median = median(dim)) )$median
  
  # computes linear regression of band PGS accuracy against the median band distance
  # uses sample size of band as weight
  lm1 <- lm(R2 ~ median, data = band_data, weights = band_data$N)
  m <- lm1$coefficients[[2]]
  m_se <- summary(lm1)$coefficients[2,2]
  
  # gets PGS accuracy of band closest to zero (by median)
  i_ref <- order(band_data$median)[1]
  R2_ref <- band_data$R2[i_ref]
  # adjust m slope by the r_ref and range
  m_hat <- m * diff(range.) / R2_ref %>% unname()
  m_hat_se <- m_se * diff(range.) / R2_ref %>% unname()
  
  # computes correlation between phenotype and PGS for entire sample
  cor1 <- cor.test(data$pheno, data$PGS)
  
  # saves data to output list
  output$band_data <- band_data
  if (!return_objects) {
    output$lm$intercept<- lm1$coefficients[[1]]
    output$lm$m <- m
    output$lm$m_se <- m_se
    output$lm$p <- summary(lm1)$coefficients[2,4]
  } else {output$lm <- lm1}
  output$lm$m_hat <- m_hat
  output$lm$m_hat_se <- m_hat_se
  if (!return_objects) {
    output$global$r <- cor1$estimate
    output$global$R2 <- cor1$estimate^2
    output$global$p <- cor1$p.value
  } else {output$global <- cor1}
  
  
  return(output)
  
}