# generates fake data represeing GxE model

library(tidyverse)

N <- 16000 # total number of individuals in test data set
L <- 32 # total number of dimensions in data set
full_data <- tibble(IID = 1:N)

# true phenotype parameters
B0 <- 0
B1 <- 2
B2 <- 2
H2 <- 0.5

# simulation parameters
AF <- 0.5
var_B1g1 <- B1^2 * AF * (1-AF)
var_B2g1E1 <- B2^2 * ( (AF^2 * 1) + (0^2 * AF*(1-AF)) + (AF * (1-AF) * 1^2) )
C0 <- sqrt((var_B1g1 + var_B2g1E1) * ((1-H2)/H2) )

#set.seed(2016)
full_data$g1 <- as.integer(sample(c(0,1),N,replace=TRUE, prob=c(1-AF,AF)))
full_data$E1 <- rnorm(N)
full_data$err <- rnorm(N)
full_data$Y <- B0 + (B1 * full_data$g1) + (B2 * full_data$g1 * full_data$E1) + (C0 * full_data$err)
#full_data$Y <- B0 + (B1 * full_data$g1) + (C0 * full_data$err)

#ggplot(full_data, aes(x=E1, y=Y, color=as.factor(g1))) + geom_point(alpha=0.05)

true_model <- lm(Y ~ g1 + g1*E1, data=full_data)
summary(true_model)


simple_model <- lm(Y ~ g1, data=full_data)
summary(simple_model)

predict_data <- full_data %>%
  mutate(Y_hat = predict(simple_model),
         SE = (Y_hat - Y)^2)

ggplot(predict_data, aes(x=abs(E1), y=SE)) +
  geom_point(alpha = 0.05) +
  geom_smooth(method='lm')

lm_error <- lm(SE ~ abs(E1),data=predict_data)
summary(lm_error)

# frNN
library(dbscan)

R <- 0.25
data_dim <- full_data %>% select(E1)
t1 <- Sys.time()
NN1 <- dbscan::frNN(data_dim, eps=R)
Sys.time() - t1

# calculates local PGS accuracy
locPGSaccs <- c()
#origin_dists <- c()
for (i in 1:N) {
  #origin_dist <- sum(data[i,2:(1+L)]**2)
  #origin_dists <- c(origin_dists, origin_dist)
  
  neighbors <- c(i,NN1$id[[i]])
  local_data <- predict_data[neighbors,c("Y_hat", "Y")]
  cor1 <- cor(local_data$Y_hat,local_data$Y)
  locPGSaccs <- c(locPGSaccs, cor1**2)
  if (i %% 100 == 0) {print(i)}
}
predict_data$locPGSacc <- locPGSaccs

ggplot(predict_data, aes(x=E1, y=locPGSaccs)) +
  geom_point(alpha = 0.05) +
  geom_smooth(method='lm') +
  ylim(0,1)

ggplot(local_data, aes(x=Y,y=Y_hat)) +
  geom_point(alpha=0.1)
