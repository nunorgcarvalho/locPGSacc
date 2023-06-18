setwd("~/group_nuno/locPGSacc/code/")
source("../code/locPGSacc.R")

if (!exists("data")) {
  data_raw <- as_tibble(fread("~/group_nuno/locPGSacc/scratch/pheno.txt"))
  N <- 10000
  i_keep <- sample(1:nrow(data_raw),N, replace = FALSE)
  data <- data_raw %>% filter(row_number() %in% i_keep) %>% select(-BMI_z)
}

col_dims <- paste0("pc",1:40)
col_pheno <- "BMI"
col_PGS <- "BMI_PGS"
R <- 25
k <- 100


data_output <- locPGSacc(data,
                         col_dims = col_dims,
                         col_pheno = col_pheno,
                         col_PGS = col_PGS,
                         #R = R,
                         k = k,
                         fixed_radius=FALSE)

colnames(data_output)
hist(log2(data_output$n_neighbors))
ggplot(data_output, aes(x=pc1, y=locPGSacc)) +
  geom_point(alpha=0.1) + geom_smooth(method='lm')
# ggplot(data_output, aes(x=n_neighbors, y=locPGSacc)) +
#   geom_point(alpha=0.1) + geom_smooth(method='lm') +xlim(200,2000)

dim_means <- colMeans(data_output %>% select(all_of(col_dims)))
data_output <- data_output %>%
  mutate(PC_dist = sqrt(rowSums((select(., all_of(col_dims)) - dim_means)^2)))
ggplot(data_output[data_output$PC_dist<50,], aes(x=PC_dist, y=locPGSacc)) +
  geom_point(alpha=0.02) + geom_smooth(method='lm')
