#' @title get_locPGS_decay
#' @description Returns statistics related to local PGS accuracy as a function of sample distance (NN method)
#' @inherit locPGSacc author
#' 
#' @details Takes a dataset with columns for local PGS accuracy (or any other
#' desired metric) and sample distance (e.g. genetic PC distance) and returns the
#' decay of PGS accuracy as a function of distance. Other optional arguments
#' allow for plotting of PGS accuracy across entire sample and within groups
#' (e.g. ancestries)
#' 
#' @param data data table containing all necessary columns and rows
#' @param col_PGSacc character: column name of the local PGS accuracy
#' @param col_dist character: column name of sample's distance in space
#' @param dist_limits (optional) numeric vector: minimum and maximum x-axis limits for distance
#' @param col_pheno (optional) character: column name of the phenotype of interest
#' @param col_PGS (optional) character: column name of the polygenic scores for the phenotype of interest
#' @param col_group (optional) character: column name of sample group assignments
#' @param col_n_neighbors (optional) character: column name of points' number of neighbors
#' @param return_objects (optional) logical: whether the function should return the outputs of cor.test() and lm() directly, rather than extracting the most important metrics
#' @param ref_window (optional) numeric: proportion of sample with lowest 'col_dist' that is used for computing standardized PGS decay slope (m_hat)
#' 
#' @return Returns a nested list with statistics related to PGS decay:
#' \itemize{
#'   \item cor: computes correlation between locPGSacc & dim_dist
#'    \itemize{
#'      \item r = Pearson's correlation coefficient
#'      \item p = p-value
#'      \item CI95 = 95% confidence interval
#'    }
#'   \item lm: computes a linear regression for locPGSacc ~ dim_dist
#'    \itemize{
#'      \item intercept = intercept of line of best fit
#'      \item m = slope of line of best fit
#'      \item p = p-value for slope
#'      \item m_hat = standardized slope; m divided by the cor(pheno, PGS) among
#'            reference population (the ref_window*100% individuals (default 5%)
#'            with the lowest 'col_dim' in the data); uses mean
#'            locPGSacc when columns not given
#'    }
#'   \item global: computes cor(pheno, PGS) on all samples
#'    \itemize{
#'      \item r = Pearson's correlation coefficient
#'      \item p = p-value
#'    }
#'   \item group: tibble containing series of statistics for each group (excluding NA) in dataset
#'    \itemize{
#'      \item group = name of group
#'      \item N = number of total individuals in group
#'      \item N_anchors = number of anchors for which locPGSacc was computed
#'      \item mean_neighbors = average number of neighbors among groups' anchors
#'      \item dist_mean = average dim_dist of all individuals in group
#'      \item dist_sd = standard deviation of dim_dist of all individuals in group
#'      \item r = correlation between phenotype and PGS for all individuals in group
#'      \item p = p-value of cor(pheno, PGS)
#'      \item r_rel = cor(pheno,PGS) of group divided by cor(pheno,PGS) among
#'            reference population (see above)
#'      \item mean_anchor_acc = average locPGSacc among groups' anchors
#'    }
#' } 
#' 
#' @export
#' 
#' @import tidyverse

get_locPGS_decay <- function(
    data,
    col_PGSacc = "locPGSacc",
    col_dist = "dim_dist",
    dist_limits = NA, # optional
    col_pheno = NA, #optional
    col_PGS = NA, #optional
    col_group = NA, # optional
    col_n_neighbors = NA,
    return_objects = FALSE,
    ref_window = 0.1 # central proportion of dim_dist range considered for standardizing m
) {
  # checks if column names exist ####
  cols <- c(col_PGSacc, col_dist, col_group, col_pheno, col_PGS, col_n_neighbors)
  cols <- cols[!is.na(cols)]
  if (!all(cols %in% colnames(data))) {stop("At least one column not found in inputted data")}
  
  # output list is established ####
  output <- list(
    cor = list(
      r = as.numeric(NA),
      p = as.numeric(NA),
      CI95 = as.numeric(NA)
    ),
    lm = list(
      intercept = as.numeric(NA),
      m = as.numeric(NA),
      p = as.numeric(NA)
    ),
    global = list(
      r = as.numeric(NA),
      p = as.numeric(NA) 
    ),
    group = tibble(group = as.character(),
                   N = as.numeric(),
                   N_anchors = as.numeric(),
                   mean_neighbors = as.numeric(),
                   dist_mean = as.numeric(),
                   dist_sd = as.numeric(),
                   r = as.numeric(),
                   p = as.numeric(),
                   r_rel = as.numeric(),
                   mean_anchor_acc = as.numeric()
                   )
  )
  # renames known columns for data set ####
  data <- data  %>% rename(locPGSacc = !!sym(col_PGSacc),
                           dim_dist = !!sym(col_dist)   )
  if (!is.na(col_pheno) & !is.na(col_PGS)) {
    data <- data %>% rename(pheno = !!sym(col_pheno),
                            PGS = !!sym(col_PGS)    )
  }
  if (!is.na(col_group)) {
    data <- data %>% rename(group = !!sym(col_group) )
  }
  if (!is.na(col_n_neighbors)) {
    data <- data %>% rename(n_neighbors = !!sym(col_n_neighbors) )
  }
  
  # filtering dataset ####  
  # gets dataset for just the anchors
  data_anchors <- data %>% filter(!is.na(locPGSacc), !is.na(dim_dist))
  # gets dist_limits and imposes them
  if (any(is.na(dist_limits))) {dist_limits <- range(data_anchors$dim_dist)}
  data_anchors <- data_anchors %>% filter(between(dim_dist, dist_limits[1], dist_limits[2]))
  
  # correlation ####
  # computes correlation statistics for locPGSacc ~ dim_dist
  cor1 <- cor.test(data_anchors$locPGSacc,
                   data_anchors$dim_dist  )
  if (!return_objects) {
    output$cor$r <- cor1$estimate
    output$cor$p <- cor1$p.value
    output$cor$CI95 <- cor1$conf.int[1:2]
  } else {output$cor <- cor1}
  
  # linear regression ####
  # computes linear regression statistics for locPGSacc ~ dim_dist
  lm1 <- lm(locPGSacc ~ dim_dist, data = data_anchors)
  m <- lm1$coefficients[[2]]
  if (!return_objects) {
    output$lm$intercept<- lm1$coefficients[[1]]
    output$lm$m <- m
    output$lm$p <- summary(lm1)$coefficients[2,4]
  } else {output$lm <- lm1}
  ## gets standardized m ####
  # range2 <- window_prop * diff(dist_limits)
  # dist_limits2 <- c(mean(dist_limits) - range2/2,
  #                   mean(dist_limits) + range2/2)
  # uses actual cor(PGS, pheno) for standard ref if given
  if (!is.na(col_pheno) & !is.na(col_PGS)) {
    # data_ref <- data %>% filter(dim_dist < dist_limits2[1],
    #                             !is.na(pheno), !is.na(PGS))
    data_ref <- data %>% filter(!is.na(pheno), !is.na(PGS))
    data_ref <- data_ref %>% arrange(dim_dist) %>%
      filter(row_number() <= ceiling(ref_window*nrow(data_ref)))
    r_ref <- cor(data_ref$pheno, data_ref$PGS)
  } else {
    #data_anchors_ref <- data_anchors %>% filter(dim_dist < dist_limits2[1])
    data_anchors_ref <- data_anchors %>% arrange(dim_dist) %>%
      filter(row_number() <= ceiling(ref_window*nrow(data_anchors)))
    r_ref <- mean(data_anchors_ref$locPGSacc)
  }
  
  m_hat <- m / r_ref
  output$lm$m_hat <- m_hat
  
  # accuracy ####
  # checks if user provided column name for phenotype and for PGS
  if (!is.na(col_pheno) & !is.na(col_PGS)) {
    ## global ####
    # computes global correlation between phenotype and PGS
    data_cor <- data %>% filter(!is.na(pheno), !is.na(PGS))
    cor2 <- cor.test(data_cor$pheno, data_cor$PGS)
    if (!return_objects) {
      output$global$r <- cor2$estimate
      output$global$p <- cor2$p.value
    } else {output$global <- cor2}
    
    # groups ####
    if (!is.na(col_group)) {
      # extracts groups and removes NA's
      groups <- data_cor$group %>% unique()
      groups <- groups[!is.na(groups)]
      for (group. in groups) {
        # computes correlation between phenotype and PGS for specific group
        data_group <- data_cor %>% filter(group == group.)
        if (nrow(data_group) <= 2) {next}
        cor3 <- cor.test(data_group$pheno, data_group$PGS)
        
        data_anchors_group <- data_anchors %>% filter(group == group.)
        #if (nrow(data_anchors_group)==0) {mean_anchor_acc <- NA
        #} else {mean_anchor_acc <- mean(data_anchors_group$locPGSacc)}
        if (is.na(col_n_neighbors)) {mean_neighbors <- NA
        } else {mean_neighbors <- mean(data_anchors_group$n_neighbors)}
        
        # adds to table
        output$group <- output$group %>%
          add_row(
            group = group.,
            N = nrow(data_group),
            N_anchors = nrow(data_anchors_group),
            mean_neighbors = mean_neighbors,
            dist_mean = mean(data_group$dim_dist),
            dist_sd = sd(data_group$dim_dist),
            r = cor3$estimate,
            p = cor3$p.value,
            r_rel = cor3$estimate / r_ref,
            mean_anchor_acc = mean(data_anchors_group$locPGSacc)
          )
      }
      # sorts groups by dist_mean
      output$group <- output$group %>%
        arrange(dist_mean)
    } else {output$group <- NULL} # removes if columns not present
    
  } else {output$global <- NULL} # removes if columns not present
  
  # returns output ####
  return(output)
}
