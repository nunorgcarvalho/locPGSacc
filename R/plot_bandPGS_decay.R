#' @title plot_bandPGS_decay
#' @description Computes PGS accuracy along bands of another variable
#' @inherit locPGSacc author
#' 
#' @details Plots the decay of PGS accuracy along some dimensional variables by
#' cutting up the dimensional variable into bands/sections. More details in
#' [get_bandPGS_decay()].
#' 
#' @inheritParams get_bandPGS_decay
#' @param fixed_ymin (optional) logical: whether the y-axis should start at 0 (default) or just be the minimum accuracy
#' 
#' @return Plot of the results of [get_bandPGS_decay()].
#' 
#' @export
#' 
#' @import tidyverse

plot_bandPGS_decay <- function (
    data,
    col_dim,
    col_pheno,
    col_PGS,
    i_omit = c(),
    window = 0.95,
    bands = 10,
    min_samples = 3,
    fixed_ymin = TRUE
) {
  
  output <- get_bandPGS_decay(data = data,
                              col_dim = col_dim,
                              col_pheno = col_pheno, # for r
                              col_PGS = col_PGS, # for r
                              i_omit = i_omit,
                              window = window,
                              bands = bands,
                              min_samples = min_samples)
  
  band_data <- output$band_data
  
  pval_text <- pvalue2text(output$lm$p)
  m_text <- paste0(round(output$lm$m*1000,3),"%*%10^-3")
  m_hat_text <- paste0(round(output$lm$m_hat*1000,3),"%*%10^-3")
  annotation <- paste0("p==",pval_text,
                       "~~hat(m)==",m_hat_text,
                       "~~m==",m_text)
  
  gg <- ggplot(band_data %>% filter(N >= 30), aes(x = median, y = r)) +
    geom_hline(yintercept = output$global$r, color="gray10") +
    geom_point() +
    geom_errorbar(aes(ymin = r_lower, ymax = r_upper)) +
    geom_abline(intercept = output$lm$intercept, slope = output$lm$m, color = "blue")
  
  # extracts x and y scale from plot so far
  xrange <- layer_scales(gg)$x$range$range
  yrange <- layer_scales(gg)$y$range$range
  
  gg <- gg +
    annotate("text", x = max(xrange), y = max(yrange), label = annotation,
             parse = TRUE, vjust=1, hjust=1) +
    #geom_smooth(method = "lm", aes(weight = N)) +
    theme_light() +
    xlab("Distance") + ylab("PGS Accuracy")
  
  if (fixed_ymin) {gg <- gg + ylim(min(0,min(yrange)), max(yrange))}
  
  return(gg)
}