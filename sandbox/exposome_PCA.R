# exploring PCs of the exposome
library(tidyverse)
library(data.table)
library(ggbiplot)

loc_pheno <- "~/scratch3/PXS_pipeline/pheno_EC.txt"
pheno <- as_tibble(fread(loc_pheno))

pheno_matrix <- pheno %>% select(starts_with("f"), -FID) %>%
  drop_na()

pca1 <- prcomp(pheno_matrix)

pca1_tibble <- as_tibble(pca1$x)

ggplot(pca1_tibble, aes(x=PC1, y=PC2)) +
  geom_point(alpha=0.05)


pca1_var <- as_tibble(t((summary(pca1))$importance)) %>%
  mutate(PC = row_number())
ggplot(pca1_var, aes(x=PC, y=`Proportion of Variance`)) +
  geom_col(width=1) +
  geom_line(mapping=aes(y=`Cumulative Proportion`*max(pca1_var$`Proportion of Variance`))) +
  scale_y_continuous(
    name = "Proportion of Variance",
    sec.axis = sec_axis( trans=~.*(1/max(pca1_var$`Proportion of Variance`)), name="Cumulative Variance")
  )

pheno_matrix2 <- pheno_matrix %>% mutate(PC1 = pca1_tibble$PC1)

lm1 <- lm(PC1 ~ ., data=pheno_matrix2)
summary(lm1)

####
# interpolating, subset
library(sinkr)

pheno_subset <- pheno %>% select(starts_with("f"), -FID) %>%
  filter(row_number() %in% sample(1:nrow(pheno),40000,replace=FALSE))

dineof1 <- dineof(as.matrix(pheno_subset))
pheno_subset_imputed <- as_tibble(dineof1$Xa)

pca2 <- prcomp(pheno_subset_imputed)

pca2_tibble <- as_tibble(pca2$x)
summary(pca2)

# interpolating, full
pheno_matrix3 <- pheno %>% select(starts_with("f"), -FID)
dineof2 <- dineof(as.matrix(pheno_matrix3))

pheno_imputed <- as_tibble(dineof2$Xa)

pca3 <- prcomp(pheno_imputed)

pca3_tibble <- as_tibble(pca3$x)
summary(pca3)
