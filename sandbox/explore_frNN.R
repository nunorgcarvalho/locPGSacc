# explore_knn - working with KNN function from exisitng package 
library(tidyverse)
library(data.table)
library(dbscan)

full_data <- as_tibble(fread("full_test_data.txt"))
col_prefix <- "PC"
seed <- 2016

N <- 8000 
L <- 32
R <- 0.25

col_dims <- paste0(col_prefix,1:L)
set.seed(seed)
data <- full_data %>%
  filter(row_number() %in% sample(1:nrow(full_data),N,replace = FALSE))

data_dim <- data %>% select(all_of(col_dims))
t1 <- Sys.time()
NN1 <- dbscan::frNN(data_dim, eps=R)
Sys.time() - t1

mean(sapply(NN1$id, length))
median(sapply(NN1$id, length))
hist(sapply(NN1$id, length))

# ggplot(data, aes(x=PC1, y=PC2)) +
#   geom_point(alpha=0.05) +
#   geom_point(data=data[NN1$id[[1]],], color="red", alpha=0.25) +
#   geom_point(data=data[1,], color="red") +
#   coord_fixed()

# calculates local PGS accuracy
locPGSaccs <- c()
#origin_dists <- c()
for (i in 1:N) {
  #origin_dist <- sum(data[i,2:(1+L)]**2)
  #origin_dists <- c(origin_dists, origin_dist)
  
  neighbors <- c(i,NN1$id[[i]])
  local_data <- data[neighbors,c("PGS", "phenotype")]
  cor1 <- cor(local_data$PGS,local_data$phenotype)
  locPGSaccs <- c(locPGSaccs, cor1**2)
  if (i %% 100 == 0) {print(i)}
}
data$locPGSacc <- locPGSaccs
#data$origin_dist <- origin_dists

# local PGS as a function of the PCs
ggplot(data, aes(x=PC1, y=locPGSaccs)) +
  geom_point(alpha = 0.1)
ggplot(data, aes(x=PC2, y=locPGSaccs)) +
  geom_point(alpha = 0.1)
ggplot(data, aes(x=PC3, y=locPGSaccs)) +
  geom_point(alpha = 0.1)
ggplot(data, aes(x=PC4, y=locPGSaccs)) +
  geom_point(alpha = 0.1)
ggplot(data, aes(x=origin_dist, y=locPGSaccs)) +
  geom_point(alpha = 0.1)

data_lm <- data %>% select(locPGSacc, starts_with(col_prefix)) %>%
  mutate_at(vars(starts_with("PC")), abs)
lm1 <- lm(locPGSacc ~ ., data=data_lm)
summary(lm1)

lm2 <- lm(locPGSacc ~ origin_dist, data=data)
summary(lm2)
