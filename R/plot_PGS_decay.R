#' @title plot_PGS_decay
#' @description Plots local PGS accuracy as a function of sample distance
#' @inherit locPGSacc author
#' 
#' @details Takes a dataset with columns for local PGS accuracy (or any other
#' desired metric) and sample distance (e.g. genetic PC distance) and plots the
#' decay of PGS accuracy as a function of distance. Other optional arguments
#' allow for plotting of PGS accuracy across entire sample and within groups
#' (e.g. ancestries). Annotated statistics are computed from [get_PGS_decay()].
#' 
#' @inheritParams get_PGS_decay
#' 
#' @export
#' 
#' @import tidyverse


plot_PGS_decay <- function(
    data,
    col_PGSacc = "locPGSacc",
    col_dist = "dim_dist",
    col_group = NA, # optional
    dist_limits = NA, # optional
    col_pheno = NA, #optional
    col_PGS = NA #optional
) {
  data_full <- data
  data <- data[!is.na(data[[col_PGSacc]]),]
  # checks if columns are in data table
  if (!col_PGSacc %in% colnames(data)) {stop(paste0("'",col_PGSacc,"' column not in data table. Use locPGSacc() or locPGSacc.FAST() functions"))}
  if (!col_dist %in% colnames(data)) {stop(paste0("'",col_dist,"' column not in data table. Use dim_dist() function."))}
  # manually sets dist_limits to full range if unspecified
  if (any(is.na(dist_limits))) {dist_limits <- range(data[,col_dist])}
  if (length(dist_limits) != 2) {stop("Length of dist_limits not equal to 2")}
  
  # makes new data table with just relevant columns and data
  data_plot <- data %>%
    dplyr::rename(dist = !!enquo(col_dist),
                  locPGSacc = !!enquo(col_PGSacc)) %>%
    filter(dist >= min(dist_limits), dist <= max(dist_limits))
  # extracts number of remaining data points and sets scatterplot alpha accordingly
  n_data_points <- sum(!is.na(data_plot$locPGSacc))
  alpha <- min(1, 2000 / n_data_points)
  
  # beginning of ggplot
  gg <- ggplot(data_plot, aes(x = dist, y = locPGSacc)) +
    theme_light()
  
  # gets decay statistics from get_PGS_decay()
  output <- get_PGS_decay(data_full,
                          col_PGSacc = col_PGSacc,
                          col_dist = col_dist,
                          dist_limits = dist_limits,
                          col_pheno = col_pheno,
                          col_PGS = col_PGS,
                          col_group = col_group)
  
  # adds horizontal lines denoting PGS accuracy within each defined group (if given)
  if (!is.na(col_group)) {
    # makes temporary ggplot for extracting resulting groups and group colors
    gg1 <- gg + geom_point(alpha = alpha, aes(color = !!sym(col_group)))
    ggb <- ggplot_build(gg1)
    group_colors <- ggb$data[[1]]$colour %>% unique()
    groups <- data_plot[[col_group]] %>% unique()
    # throws error if columns for the phenotype and PGS were not given
    if (is.na(col_pheno) | is.na(col_PGS)) {stop("'col_group' was given but 'col_pheno' and 'col_PGS' were not")}
    
    for (i in 1:length(groups)) {
      group. <- groups[i]
      if (is.na(group.)) {next}
      group_color <- group_colors[i]
      # gets mean of locPGSacc among group's anchors
      mean_anchor_acc <- (output$group %>%
        filter(group == group.) )$mean_anchor_acc
      # adds horizontal lines corresponding to r in each group
      gg <- gg + geom_hline(yintercept = mean_anchor_acc, color=group_color)
    }
    # adds points of samples, colored by their group
    gg <- gg + geom_point(alpha = alpha, aes(color = !!sym(col_group)))
  } else {
    # adds points of samples
    gg <- gg + geom_point(alpha = alpha)
  }
  
  # if columns for phenotype and PGS were given, computes cor(phenotype, PGS) for entire sample
  if (!is.na(col_pheno) & !is.na(col_PGS)) {
    r_global <- output$global$r
    gg <- gg + geom_hline(yintercept = r_global, color="grey")
  }
  
  # gets statistics to annotate 
  r <- output$cor$r
  p <- output$cor$p
  pval_text <- pvalue2text(p)
  m <- output$lm$m
  m_text <- paste0(round(m*1000,3),"%*%10^-3")
  m_hat <- output$lm$m_hat
  m_hat_text <- paste0(round(m_hat*1000,3),"%*%10^-3")
  # pre-writes annotation of locPGS ~ dist using plotmath() parsing
  annotation <- paste0("r==",round(r,3),"~~p==",pval_text,"~~hat(m)==",m_hat_text,"~~m==",m_text)
  
  # extracts x and y scale from plot so far
  xrange <- layer_scales(gg)$x$range$range
  yrange <- layer_scales(gg)$y$range$range
  
  # makes rest of points
  gg <- gg + 
    geom_smooth(method='lm', color="red", formula = y ~ x) +
    annotate("text", x = max(xrange), y = max(yrange), label = annotation,
             parse = TRUE, vjust=1, hjust=1) +
    xlab("Distance") + ylab("Local PGS Accuracy") +
    labs(color="Group") +
    guides(color = guide_legend(override.aes = list(alpha = 1))) # makes legend color alpha = 1
  
  # returns ggplot
  return(gg)
}