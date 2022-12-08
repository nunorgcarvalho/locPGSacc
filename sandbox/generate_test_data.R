# generates fake PC-like data to test on

library(tidyverse)

N <- 16000 # total number of individuals in test data set
L <- 32 # total number of dimensions in data set


set.seed(2016)
full_data <- tibble(IID = 1:N)
for (j in 1:L) {
  # data for dimensions is normally distributed
  # dimension j's variance is 1/j
  full_data[,paste0("PC",j)] <- rnorm(N, mean=0, sd = (1/j)**2)
  print(j)
}

# generates fake output data based on PC data
PC_coeffs <- rnorm(L)
full_data$PGS <- as.numeric(NA)
full_data$phenotype <- rnorm(N)
PC_error_prop <- 0.5
absolute_error_prop <- 1
for (i in 1:nrow(data)) {
  pos <- data[i,2:(1+L)] %>% unlist() %>% unname()
  PC_error <- PC_coeffs %*% pos
  #origin_dist <- dist(rbind(data[i,2:(1+L)], rep(0,L)))[[1]]
  PGS <- full_data$phenotype[i] + PC_error_prop*(PC_error) + absolute_error_prop*rnorm(1)
  full_data$PGS[i] <- PGS
  if (i %% 100 == 0) {print(i)}
}
# plots (PGS-phenotype) by PC1
ggplot(full_data, aes(x=PC3, y=(PGS-phenotype))) +
  geom_point(alpha=0.03) +
  geom_smooth(method="lm", color="blue") +
  geom_abline(slope=PC_coeffs[3], intercept = 0)

# plots PGS vs phenotype
ggplot(full_data, aes(x=phenotype, y=PGS)) +
  geom_point(alpha=0.05)
cor.test(full_data$PGS, full_data$phenotype)

# plots PC data for first two dimensions
ggplot(full_data, aes(x=PC1, y=PC2)) +
  geom_point(alpha = 0.1)

write.table(full_data, "full_test_data.txt", sep="\t",row.names=FALSE, quote=FALSE)
