# generates fake PC-like data to test on

library(tidyverse)

N <- 16000 # total number of individuals in test data set
L <- 32 # total number of dimensions in data set


set.seed(2016)
data <- tibble(IID = 1:N)
for (i in 1:L) {
  # data for dimensions is normally distributed
  # dimension i's variance is 1/i
  data[,paste0("PC",i)] <- rnorm(N, mean=0, sd = (1/i)**2)
  print(i)
}

# plots data for first two dimensions
ggplot(data, aes(x=PC1, y=PC2)) +
  geom_point(alpha = 0.1)

write.table(data, "full_test_data.txt", sep="\t",row.names=FALSE, quote=FALSE)
