# plots locPGSacc against PC distance
library(tidyverse)
source("R/pvalue_to_text.R")


plot_PGS_decay <- function(
    data,
    col_PGSacc = "locPGSacc",
    col_dist = "dim_dist",
    #col_n_neighbors = "n_neighbors", # set to "" if not wanted to be plotted
    col_group = NA, # optional
    dist_limits = NA, # optional
    col_pheno = NA, #optional
    col_PGS = NA #optional
    
) {
  # checks if columns are in data table
  if (!col_PGSacc %in% colnames(data)) {stop(paste0("'",col_PGSacc,"' column not in data table"))}
  if (!col_dist %in% colnames(data)) {stop(paste0("'",col_dist,"' column not in data table"))}
  # manually sets dist_limits to full range if unspecified
  if (any(is.na(dist_limits))) {dist_limits <- range(data[,col_dist])}
  if (length(dist_limits) != 2) {stop("Length of dist_limits not equal to 2")}
  
  # handles n_neighbors column depending on input
  # if (col_n_neighbors == "") {
  #   col_n_neighbors <- "n_neighbors"
  #   data[,col_n_neighbors] <- -1
  # } else {
  #   if (!col_n_neighbors %in% colnames(data)) {stop(paste0("'",col_n_neighbors,"' column not in data table"))}
  # }
  
  # makes new data table with just relevant columns and data
  data_plot <- data %>%
    dplyr::rename(dist = !!enquo(col_dist),
                  locPGSacc = !!enquo(col_PGSacc)) %>%
    filter(dist >= min(dist_limits), dist <= max(dist_limits))
  n_data_points <- sum(!is.na(data_plot$locPGSacc))
  alpha <- min(1, 2000 / n_data_points)
  
  gg <- ggplot(data_plot, aes(x = dist, y = locPGSacc)) +
    theme_light()
  
  if (!is.na(col_group)) {
    gg1 <- gg + geom_point(alpha = alpha, aes(color = !!sym(col_group)))
    ggb <- ggplot_build(gg1)
    group_colors <- ggb$data[[1]]$colour %>% unique()
    #groups_i <- ggb$data[[1]]$group %>% unique()
    groups <- data_plot[[col_group]] %>% unique()
    
    if (!is.na(col_pheno) & !is.na(col_PGS)) {
      
      for (i in 1:length(groups)) {
        group <- groups[i]
        if (is.na(group)) {next}
        group_color <- group_colors[i]
        # r_group <- cor(data_plot[which(data_plot[col_group]==group & !is.na(data_plot[col_pheno])),][[col_pheno]],
        #                 data_plot[which(data_plot[col_group]==group & !is.na(data_plot[col_PGS])),][[col_PGS]])
        r_group <- mean(data_plot[data_plot[[col_group]]==group,][[col_PGSacc]], na.rm = TRUE)
        gg <- gg + geom_hline(yintercept = r_group, color=group_color)
      }
    }
    gg <- gg + geom_point(alpha = alpha, aes(color = !!sym(col_group)))
  } else {
    gg <- gg + geom_point(alpha = alpha)
  }
  
  if (!is.na(col_pheno) & !is.na(col_PGS)) {
    r_global <- cor(data_plot[[col_pheno]], data_plot[[col_PGS]])
    gg <- gg + geom_hline(yintercept = r_global, color="grey")
  }
  
  cor1 <- cor.test(data_plot$locPGSacc, data_plot$dist)
  lm1 <- lm(locPGSacc ~ dist, data = data_plot)
  
  r <- round(cor1$estimate,3)
  pval_text <- pvalue_to_text(cor1$p.value)
  slope <- pvalue_to_text(lm1$coefficients[[2]])
  
  
  xrange <- layer_scales(gg)$x$range$range
  yrange <- layer_scales(gg)$y$range$range
  annotation <- paste0("r==",r,"~'\n'~","~p==",pval_text,"~'\n'~","m==",slope)
  gg <- gg + 
    geom_smooth(method='lm', color="red") +
    annotate(
      "text", x = max(xrange), y = max(yrange),
      label = annotation, parse = TRUE,
      vjust=1, hjust=1
      ) +
    xlab("Distance") +
    ylab("Local PGS Accuracy") +
    labs(color="Group") +
    guides(color = guide_legend(override.aes = list(alpha = 1)))
  gg
}