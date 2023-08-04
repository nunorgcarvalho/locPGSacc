#' @title get_bandPGS_decay
#' @description Returns statistics related to local PGS accuracy as a function of sample distance (band method)
#' @inherit locPGSacc author
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

get_bandPGS_decay <- function (
    data,
    col_dim,
    col_pheno, # for r
    col_PGS, # for r
    i_omit = c(),
    window = 0.95,
    bands = 10
) {
  
  # output list is established ####
  output <- list(
    band_data = tibble(),
    lm = list(
      intercept = as.numeric(NA),
      m = as.numeric(NA),
      p = as.numeric(NA)
    ),
    global = list(
      r = as.numeric(NA),
      p = as.numeric(NA) 
    )
  )
  
  
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
  
  band_data <- tibble(band = as.character(),
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
    
    band_data <- band_data %>% add_row(
      band = band,
      min = min(data_band$dim),
      max = max(data_band$dim),
      N = nrow(data_band),
      r = cor1$estimate,
      r_lower = cor1$conf.int[1],
      r_upper = cor1$conf.int[2]
    )
  }
  band_data$median <- ( data %>% group_by(dim_group) %>% summarize(median = median(dim)) )$median
  lm1 <- lm(r ~ median, data = band_data, weights = band_data$N)
  cor1 <- cor.test(data$pheno, data$PGS)
  
  
  output$band_data <- band_data
  output$lm$intercept<- lm1$coefficients[[1]]
  output$lm$m <- lm1$coefficients[[2]]
  output$lm$p <- summary(lm1)$coefficients[2,4]
  
  if (!is.na(col_pheno) & !is.na(col_PGS)) {
    output$global$r <- cor1$estimate
    output$global$p <- cor1$p.value
  } else {output$global <- NULL}
  
  return(output)
  
}