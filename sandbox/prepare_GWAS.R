library(tidyverse)
library(data.table)

loc_phenoEC <- "~/scratch3/PXS_pipeline/pheno_EC.txt"
pheno_full <- as_tibble(fread(loc_phenoEC))

pheno <- pheno_full %>% select(FID,IID,sex,age, assessment_center, starts_with("pc"),BMI = f21001)
loc_out <- "~/scratch3/locPGSacc/pheno_BMI.txt"
write.table(pheno, loc_out, sep="\t", quote=FALSE, row.names = FALSE)
