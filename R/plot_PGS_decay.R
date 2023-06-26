#' @title plot_PGS_decay
#' @description Returns plot of local PGS accuracy as a function of sample distance
#' @inherit locPGSacc author
#' 
#' @details Takes a dataset with columns for local PGS accuracy (or any other
#' desired metric) and sample distance (e.g. genetic PC distance) and plots the
#' decay of PGS accuracy as a function of distance. Other optional arguments
#' allow for plotting of PGS accuracy across entire sample and within groups
#' (e.g. ancestries)
#' 
#' @param data data table containing all necessary columns and rows
#' @param col_PGSacc character: column name of the local PGS accuracy
#' @param col_dist character: column name of sample's distance in space
#' @param col_group (optional) character: column name of sample group assignments
#' @param dist_limits (optional) numeric vector: minimum and maximum x-axis limits for distance
#' @param col_pheno (optional) character: column name of the phenotype of interest
#' @param col_PGS (optional) character: column name of the polygenic scores for the phenotype of interest
#' @param r_group_mode (optional) character: mode for computing group PGS accuracy values
#' \itemize{
#'   \item 'r' = calculates cor(phenotype, PGS) for all of that group's samples
#'   \item 'mean' = calculates mean PGS accuracy for that group's anchor points
#' }
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
    col_PGS = NA, #optional
    r_group_mode = "r" # 'mean' or 'r'
) {
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
  
  # adds horizontal lines denoting PGS accuracy within each defined group (if given)
  if (!is.na(col_group)) {
    # makes temporary ggplot for extracting resulting groups and group colors
    gg1 <- gg + geom_point(alpha = alpha, aes(color = !!sym(col_group)))
    ggb <- ggplot_build(gg1)
    group_colors <- ggb$data[[1]]$colour %>% unique()
    groups <- data_plot[[col_group]] %>% unique()
    # throws error if columns for the phenotype and PGS were not given
    if (!is.na(col_pheno) & !is.na(col_PGS)) {
      # loops through each group in dataset
      for (i in 1:length(groups)) {
        group <- groups[i]
        group_color <- group_colors[i]
        # skips if group is NA
        if (is.na(group)) {next}
        # computes r in each group in different ways
        if (r_group_mode == "r") {
          # calculates cor(phenotype, PGS) for all of that group's samples
          r_group <- cor(data_plot[which(data_plot[col_group]==group & !is.na(data_plot[col_pheno])),][[col_pheno]],
                          data_plot[which(data_plot[col_group]==group & !is.na(data_plot[col_PGS])),][[col_PGS]])
        } else if (r_group_mode == "mean") {
          # calculates mean PGS accuracy for that group's anchor points 
          r_group <- mean(data_plot[data_plot[[col_group]]==group,][[col_PGSacc]], na.rm = TRUE)
        } else {stop("'r_group_mode' must be one of 'r' or 'mean'")}
        # adds horizontal lines corresponding to r in each group
        gg <- gg + geom_hline(yintercept = r_group, color=group_color)
      }
    } else {stop("'col_group' was given but 'col_pheno' and 'col_PGS' were not")}
    # adds points of samples, colored by their group
    gg <- gg + geom_point(alpha = alpha, aes(color = !!sym(col_group)))
  } else {
    # adds points of samples
    gg <- gg + geom_point(alpha = alpha)
  }
  
  # if columns for phenotype and PGS were given, computes cor(phenotype, PGS) for entire sample
  if (!is.na(col_pheno) & !is.na(col_PGS)) {
    r_global <- cor(data_plot[[col_pheno]], data_plot[[col_PGS]])
    gg <- gg + geom_hline(yintercept = r_global, color="grey")
  }
  
  # gets correlation and linear regression information on locPGSacc ~ dist
  cor1 <- cor.test(data_plot$locPGSacc, data_plot$dist)
  lm1 <- lm(locPGSacc ~ dist, data = data_plot)
  r <- round(cor1$estimate,3)
  pval_text <- pvalue2text(cor1$p.value)
  slope <- pvalue_to_text(lm1$coefficients[[2]])
  
  # extracts x and y scale from plot so far
  xrange <- layer_scales(gg)$x$range$range
  yrange <- layer_scales(gg)$y$range$range
  # pre-writes annotation of locPGS ~ dist using plotmath() parsing
  annotation <- paste0("r==",r,"~'\n'~","~p==",pval_text,"~'\n'~","m==",slope)
  
  # makes rest of points
  gg <- gg + 
    geom_smooth(method='lm', color="red") +
    annotate("text", x = max(xrange), y = max(yrange), label = annotation,
             parse = TRUE, vjust=1, hjust=1) +
    xlab("Distance") + ylab("Local PGS Accuracy") +
    labs(color="Group") +
    guides(color = guide_legend(override.aes = list(alpha = 1))) # makes legend color alpha = 1
  
  # returns ggplot
  return(gg)
}