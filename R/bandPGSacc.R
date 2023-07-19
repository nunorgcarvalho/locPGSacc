#' @title bandPGSacc
#' @description Computes PGS accuracy along bands of another variable
#' @author Nuno R. G. Carvalho: \email{nunocarvalho@@gatech.edu}
#' 
#' @details WIP
#' 
#' @inheritParams get_NNs
#' @inheritParams get_accuracy
#' @param col_dim character: column name of the dimension variable for which to compute bands from
#' @param window (optional) numeric: proportion of the distribution, centered at the median, to include when making bands
#' @param bands (optional) integer: number of bands (i.e. sections) to split up the sample by according to the dimensional variable
#' 
#' @return WIP
#' 
#'
#' @export
#' 
#' @import tidyverse

bandPGSacc <- function (
    data,
    col_dim,
    col_pheno, # for r
    col_PGS, # for r
    i_omit = c(),
    window = 0.95,
    bands = 10
) {
  
  i_omit <- (1:nrow(data))[is.na(data[[col_pheno]]) | is.na(data[[col_PGS]])]
  
  data <- data[-i_omit,] %>% select(dim = !!sym(col_dim),
                                   pheno = !!sym(col_pheno),
                                   PGS = !!sym(col_PGS))
  
  low_bound_i <- round((0.5 - window/2)*length(data$dim))
  upp_bound_i <- round((0.5 + window/2)*length(data$dim))
  range. = sort(data$dim)[low_bound_i:upp_bound_i] %>% range()
  
  breaks <- c( min(data$dim), seq(range.[1], range.[2], diff(range.)/bands), max(data$dim) )
  
  data$dim_group <- cut(data$dim, breaks = breaks, include.lowest = TRUE)
  table(data$dim_group)
  
  band_decay <- tibble(band = as.character(),
                       min = as.numeric(),
                       max = as.numeric(),
                       N = as.numeric(),
                       r = as.numeric(),
                       r_lower = as.numeric(),
                       r_upper = as.numeric()
  )
  for (band in levels(data$dim_group)) {
    #print(band)
    data_band <- data %>% filter(dim_group == band)
    
    if (nrow(data_band) <= 2) {next}
    
    cor1 <- cor.test(data_band$PGS,data_band$pheno)
    
    band_decay <- band_decay %>% add_row(
      band = band,
      min = min(data_band$dim),
      max = max(data_band$dim),
      N = nrow(data_band),
      r = cor1$estimate,
      r_lower = cor1$conf.int[1],
      r_upper = cor1$conf.int[2]
    )
  }
  band_decay$median <- ( data %>% group_by(dim_group) %>% summarize(median = median(dim)) )$median
  lm1 <- lm(r ~ median, data = band_decay, weights = band_decay$N)
  summary(lm1)
  p <- summary(lm1)$coefficients[2,4]
  r_global <- cor(data$pheno, data$PGS)
  
  print(paste(col_dim, code, p < 0.05, p))
  
  ggplot(band_decay %>% filter(N >= 30), aes(x = median, y = r)) +
    geom_hline(yintercept = r_global, color="red") +
    geom_point() +
    geom_errorbar(aes(ymin = r_lower, ymax = r_upper)) +
    #geom_abline(intercept = lm1$coefficients[1], slope = lm1$coefficients[2], color = "red") +
    geom_smooth(method = "lm", aes(weight = N))
}