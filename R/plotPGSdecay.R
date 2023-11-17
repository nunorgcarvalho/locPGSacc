#' @title plotPGSdecay
#' @description Computes PGS accuracy along bins of a continuous variable
#' @inherit locPGSacc author
#' 
#' @details Plots the decay of PGS accuracy along some dimensional variables by
#' cutting up the dimensional variable into bins. More details in [PGSdecay()].
#' 
#' @inheritParams PGSdecay
#' @param min_samples (optional) integer: minimum number of samples needed in a bin to compute PGS accuracy. Must be > 2
#' @param fixed_ymin (optional) logical: whether the y-axis should start at 0 (default) or just be the minimum accuracy
#' @param show_stats (optional) logical: whether the decay slope, standardized decay slope, and corresponding p-value are shown on the plot (top right)
#' @param plot_rel_dist (optional) logical: whether the x-axis of the plot should be scaled such that the maximum distance is 1. Does not affect portability statistics.
#' 
#' @return Plot of the results of [PGSdecay()].
#' 
#' @export
#' 
#' @import tidyverse

plotPGSdecay <- function (
    data,
    col_dim,
    col_pheno,
    col_PGS,
    i_omit = c(),
    ref_window = 0.95,
    bins = 15,
    min_samples = 30,
    fixed_ymin = TRUE,
    show_stats = TRUE,
    plot_rel_dist = FALSE
) {
  
  output <- PGSdecay(data = data,
                     col_dim = col_dim,
                     col_pheno = col_pheno,
                     col_PGS = col_PGS,
                     i_omit = i_omit,
                     ref_window = ref_window,
                     bins = bins)
  
  plot_m <- output$lm$m
  bin_data <- output$bin_data
  
  if (plot_rel_dist) {
    range95 <- diff(c(min(bin_data$min),max(bin_data$max)))
    bin_data$median <- bin_data$median / range95
    plot_m <- plot_m * range95
  }
  
  gg <- ggplot(bin_data %>% filter(N >= min_samples), aes(x = median, y = R2)) +
    geom_hline(yintercept = output$global$R2, color="gray10") +
    geom_point() +
    geom_errorbar(aes(ymin = R2_lower, ymax = R2_upper)) +
    geom_abline(intercept = output$lm$intercept, slope = plot_m, color = "blue") +
    theme_light() +
    xlab("Distance") + ylab("PGS Accuracy")
  
  # extracts x and y scale from plot so far
  xrange <- layer_scales(gg)$x$range$range
  yrange <- layer_scales(gg)$y$range$range
  
  # shows slopes and p-value
  if (show_stats) {
    pval_text <- pvalue2text(output$lm$p)
    m_text <- paste0(round(output$lm$m*1000,3),"%*%10^-3")
    m_hat_text <- paste0(round(output$lm$m_hat,3))
    annotation <- paste0("p==",pval_text,
                         "~~hat(m)==",m_hat_text,
                         "~~m==",m_text)
    
    gg <- gg + annotate("text", x = max(xrange), y = max(yrange), label = annotation,
                        parse = TRUE, vjust=1, hjust=1)
    
  }
  
  if (plot_rel_dist) {gg <- gg + scale_x_continuous(limits = c(0,1), expand=c(0,0.025))}
  if (fixed_ymin) {gg <- gg + ylim(min(0,min(yrange)), max(yrange))}
  
  return(gg)
}