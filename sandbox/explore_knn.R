# explore_knn - working with KNN function from exisitng package 
library(tidyverse)
library(data.table)
library(dbscan)

full_data <- as_tibble(fread("full_test_data.txt"))
col_prefix <- "PC"
seed <- 2016

N <- 1000 
L <- 32
R <- 0.5

col_dims <- paste0(col_prefix,1:L)
set.seed(seed)
data <- full_data %>% select(all_of(col_dims)) %>%
  filter(row_number() %in% sample(1:nrow(full_data),N,replace = FALSE))

t1 <- Sys.time()
NN1 <- dbscan::frNN(data, eps=R)
Sys.time() - t1

neighbors <- NN1$id
mean(sapply(neighbors, length))
median(sapply(neighbors, length))
hist(sapply(neighbors, length))

ggplot(data, aes(x=PC1, y=PC2)) +
  geom_point(alpha=0.05) +
  geom_point(data=data[neighbors[[1]],], color="red", alpha=0.25) +
  geom_point(data=data[1,], color="red") +
  coord_fixed()
