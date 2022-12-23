library(tidyverse)
library(data.table)

loc_phenoEC <- "~/scratch3/PXS_pipeline/pheno_EC.txt"
pheno_full <- as_tibble(fread(loc_phenoEC))

pheno <- pheno_full %>%
  select(FID,IID,sex,age, assessment_center, starts_with("pc"),BMI = f21001) %>%
  mutate(group = sample(c("A","B"), nrow(pheno_full), replace=TRUE, prob = c(0.9,0.1)))
                        
loc_out <- "~/scratch3/locPGSacc/pheno_BMI_full.txt"
write.table(pheno, loc_out, sep="\t", quote=FALSE, row.names = FALSE)

pheno_train <- pheno %>% filter(group=="A")
pheno_test <- pheno %>% filter(group=="B")

loc_out <- "~/scratch3/locPGSacc/pheno_BMI_train.txt"
write.table(pheno_train, loc_out, sep="\t", quote=FALSE, row.names = FALSE)

loc_out <- "~/scratch3/locPGSacc/pheno_BMI_test.txt"
write.table(pheno_test, loc_out, sep="\t", quote=FALSE, row.names = FALSE)