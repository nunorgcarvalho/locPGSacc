library(tidyverse)
library(data.table)
library(callr)

# loads paths, directories, and definitions
source("paths_and_definitions.R")

#### Ancestry Assortment ####

# uses population-defining method from Prive et al. 2022 (and their Github code)
PCs <- as_tibble(fread(loc_40PCs))
PC_UKBB <- PCs[,paste0("pc",1:16)]
all_centers <- read.csv(
  "https://raw.githubusercontent.com/privefl/UKBB-PGS/main/pop_centers.csv",
  stringsAsFactors = FALSE)
all_sq_dist <- apply(all_centers[-1], 1, function(one_center) {
  rowSums(sweep(PC_UKBB, 2, one_center, '-')^2)
})
thr_sq_dist <- max(dist(all_centers[-1])^2) * 0.002 / 0.16
group <- apply(all_sq_dist, 1, function(x) {
  grp <- NA
  ind <- which.min(x)
  if (isTRUE(x[ind] < thr_sq_dist)) {
    grp <- all_centers$Ancestry[ind]
    # We used a more stringent cutoff for the Ashkenazi group
    if (grp == "Ashkenazi" && x[ind] > 12.5^2) grp <- NA
  }
  grp
})
table(group, exclude = NULL)
# results almost match exactly with Prive et al. 2022

# writes out file containing individuals' ancestral group
ancestry_table <- tibble(FID = PCs$FID, IID = PCs$IID, ancestry = str_replace(group, "United Kingdom","UK"))
loc_out <- "../scratch/ancestry_table.txt"
write.table(ancestry_table, loc_out, sep="\t", quote=FALSE, row.names = FALSE)

#### Phenotype extraction ####

# list of phenotype columns to keep
cols_keep1 <- c("f.eid", paste0("f.",c(31,34,54,50,21001,2443),".0.0"))
#cols_keep2 <- c("eid", paste0(c(130708),"-0.0"))

pheno1 <- as_tibble(fread(loc_pheno_full1, select = cols_keep1))
#pheno2 <- as_tibble(fread(loc_pheno_full2, select = cols_keep2))

pheno <- pheno1 %>% #left_join(pheno2, by=c("f.eid"="eid")) %>%
  select("FID"="f.eid","IID"="f.eid","sex"="f.31.0.0","age"="f.34.0.0",
         "assessment_center"="f.54.0.0","height"="f.50.0.0", "BMI"="f.21001.0.0",
         "T2D"="f.2443.0.0") %>%
  left_join(ancestry_table, by=c("FID","IID"))
pheno <- pheno %>%
  filter(!(rowSums(pheno[,sapply(pheno, class) %in% c("integer","numeric")] < 0) %in% c(1,NA))) %>%
  left_join(PCs) %>% drop_na()

# sets all non-UK to testing, and 10,000 UK to testing (rest are training)
pheno$group <- "testing"
pheno[sample((1:nrow(pheno))[pheno$ancestry=="UK"],sum(pheno$ancestry=="UK") - 10000, replace=FALSE),"group"] <- "training"
# use to check training/testing sizes by ancestry
pheno %>% group_by(ancestry, group) %>% summarize(n=n())

# saves all phenotypes, training FIDs/IIDs, and testing FIDs/IIDs
write.table(pheno, "../scratch/pheno.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(pheno %>% filter(group=="training") %>% select(FID,IID),
            "../scratch/training_FID_IID.txt",
            sep="\t", quote= FALSE, row.names = FALSE, col.names = FALSE)
write.table(pheno %>% filter(group=="testing") %>% select(FID,IID),
            "../scratch/testing_FID_IID.txt",
            sep="\t", quote= FALSE, row.names = FALSE, col.names = FALSE)
