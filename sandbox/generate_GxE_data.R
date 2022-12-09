# generates fake data represeing GxE model

library(tidyverse)

N <- 16000 # total number of individuals in test data set
L <- 32 # total number of dimensions in data set
M <- 5 # number of SNPs
Ei <- 1 # SNP number that has GxE
full_data <- tibble(IID = 1:N)

# true phenotype parameters
H2 <- 0.6 # how much variance is explained by genetics
prop_h2gxe <- 0.1 # how much of the h2 is explained by GxE (as opposed to additive effects)
B0 <- 0
#B1 <- 2
B1 <- rnorm(M)
#B2 <- 1

# simulation parameters
#AF <- 0.5
AF <- runif(M, 0.05,0.95)
var_additive <- sum(2 * B1^2 * AF * (1-AF))
var_GxE_base <- ((2*AF[Ei])^2 * 1) + (0^2 * 2 * AF[Ei]*(1-AF[Ei])) + (2 * AF[Ei] * (1-AF[Ei]) * 1^2)
B2 <- sqrt( (var_additive/var_GxE_base) * (prop_h2gxe/(1-prop_h2gxe)))
C0 <- sqrt((var_additive + var_GxE) * ((1-H2)/H2) )

#set.seed(2016)
for (i in 1:M) {
  col_g <- paste0("g",i)
  full_data[,col_g] <- sample(c(0,1,2),N,replace=TRUE, prob=c((1-AF[i])^2,2*AF[i]*(1-AF[i]),AF[i]^2))
}
#full_data$g1 <- as.integer(sample(c(0,1),N,replace=TRUE, prob=c(1-AF,AF)))
full_data$E1 <- rnorm(N)
full_data$err <- rnorm(N)
#full_data$additive_effect <- rowSums(full_data[,paste0("g",1:M)] * B1)
full_data$additive_effect <- as.matrix(full_data[,paste0("g",1:M)]) %*% B1
full_data$Y <- B0 + (full_data$additive_effect) + (B2 * full_data$g1 * full_data$E1) + (C0 * full_data$err)
#full_data$Y <- B0 + (B1 * full_data$g1) + (C0 * full_data$err)

#ggplot(full_data, aes(x=E1, y=Y, color=as.factor(g1))) + geom_point(alpha=0.05)

true_model <- lm(Y ~ g1 + g2 + g3 + g4 + g5 + g1*E1, data=full_data)
summary(true_model)


simple_model <- lm(Y ~ g1 + g2 + g3 + g4 + g5, data=full_data)
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

ggplot(predict_data, aes(x=E1, y=locPGSacc)) +
  geom_point(alpha = 0.05) +
  geom_smooth(data=predict_data%>%filter(E1 <= 0),method='lm', color="red") +
  geom_smooth(data=predict_data%>%filter(E1 >= 0),method='lm', color="blue") +
  ylim(0,1)

ggplot(predict_data, aes(x=abs(E1), y=locPGSacc)) +
  geom_point(alpha = 0.05) +
  geom_smooth(method='lm') +
  ylim(0,1)

lm_acc <- lm(locPGSacc ~ abs(E1),data=predict_data)
summary(lm_acc)

lm_acc1 <- lm(locPGSacc ~ E1,data=predict_data%>%filter(E1 <= 0))
summary(lm_acc1)
lm_acc2 <- lm(locPGSacc ~ E1,data=predict_data%>%filter(E1 >= 0))
summary(lm_acc2)

