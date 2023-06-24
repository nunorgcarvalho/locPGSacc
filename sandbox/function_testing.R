setwd("~/group_nuno/locPGSacc/code/")

if (!exists("data_raw")) {
  data_raw <- as_tibble(fread("~/group_nuno/locPGSacc/scratch/pheno.txt"))
}

N <- 100000
#N <- 484198
i_keep <- sample(1:nrow(data_raw),N, replace = FALSE)
data <- data_raw %>% filter(row_number() %in% i_keep) %>% select(-BMI_z)

col_dims <- paste0("pc",1:40)
col_pheno <- "BMI"
col_PGS <- "BMI_PGS"
R <- 40
k <- 500

t1 <- Sys.time()
source("../code/locPGSacc.R")
data_output <- locPGSacc(data,
                         col_dims = col_dims,
                         col_pheno = col_pheno,
                         col_PGS = col_PGS,
                         R = R,
                         k = k,
                         mode="hybrid"
                         )
Sys.time() - t1

hist(log2(data_output$n_neighbors))

source("../code/dim_dist.R")
data_output <- dim_dist(data_output,
                        col_dims = col_dims,
                        reference_point = 0,
                        col_dist = "PC_dist"
                        )

source("../code/plot_PGS_decay.R")
plot_PGS_decay(data_output[data_output$PC_dist<100,],
               col_dist = "PC_dist")



#### FAST function
ggplot(data_dims, aes(x=pc1, y=pc2)) +
  geom_point(alpha=0.1) +
  geom_point(data=data_dims[housing >= coverage,], color="blue") +
  geom_point(data=data_dims[anchor,], color="red")

ggplot(data_dims, aes(x=pc1, y=pc2)) +
  geom_point(alpha=1, shape=1, aes(color=housing))# +
  #geom_point(data=data_dims[housing >= coverage,], color="blue") +
  #geom_point(data=data_dims[anchor,], color="red")

ggplot(data_dims, aes(x=pc1, y=pc2)) +
  geom_point(alpha=0.1) +
  geom_point(data=data_dims[NN_ids[[1]], ], color="blue") +
  geom_point(data=data_dims[anchor,], color="red")

ggplot(data_dims, aes(x=pc1, y=pc2)) +
  geom_point(alpha=0.1) +
  geom_point(data=data_dims[homeless, ], color="blue")
