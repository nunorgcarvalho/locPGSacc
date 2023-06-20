setwd("~/group_nuno/locPGSacc/code/")

if (!exists("data_raw")) {
  data_raw <- as_tibble(fread("~/group_nuno/locPGSacc/scratch/pheno.txt"))
  N <- 10000
  i_keep <- sample(1:nrow(data_raw),N, replace = FALSE)
  data <- data_raw %>% filter(row_number() %in% i_keep) %>% select(-BMI_z)
}

col_dims <- paste0("pc",1:40)
col_pheno <- "BMI"
col_PGS <- "BMI_PGS"
R <- 40
k <- 500

source("../code/locPGSacc.R")
data_output <- locPGSacc(data,
                         col_dims = col_dims,
                         col_pheno = col_pheno,
                         col_PGS = col_PGS,
                         R = R,
                         k = k,
                         mode="hybrid"
                         )

hist(log2(data_output$n_neighbors))

source("../code/dim_dist.R")
data_output <- dim_dist(data_output,
                        col_dims = col_dims,
                        reference_point = 0,
                        col_dist = "PC_dist"
                        )
ggplot(data_output[data_output$PC_dist<100,],
       aes(x=PC_dist, y=locPGSacc)) +
  geom_point(alpha=0.1, aes(color=log2(n_neighbors))) + geom_smooth(method='lm')

source("../code/plot_PGS_decay.R")
plot_PGS_decay(data_output[data_output$PC_dist<1000,],
               col_dist = "PC_dist")
