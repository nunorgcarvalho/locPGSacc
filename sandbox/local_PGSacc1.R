library(tidyverse)
library(data.table)
library(dbscan)
library(sinkr)

dir_repo <- "~/group_nuno/locPGSacc/"
dir_scratch <- paste0(dir_repo, "scratch/")
# BMI-PGS
cols_to_keep <- c("eid",paste0(c(21001,26216),"-0.0"))
loc_pheno1 <- "/n/groups/patel/uk_biobank/project_22881_669542/ukb669542.csv"
#cols_to_keep1 <- c("eid",paste0(26216,"-0.0"))
pheno1 <- as_tibble(fread(loc_pheno1, select=cols_to_keep))
# BMI
loc_pheno2 <- "/n/groups/patel/uk_biobank/main_data_34521/ukb34521.csv"
#cols_to_keep2 <- c("eid",paste0(21001,"-0.0"))
cols_to_keep2 <- c("eid",cols_to_keep[!cols_to_keep %in% colnames(pheno1)])
pheno2 <- as_tibble(fread(loc_pheno2, select=cols_to_keep2))
# PCs
loc_PCs <- "~/scratch3/key_data/UKB_40PCs_500k.txt"
PCs <- as_tibble(fread(loc_PCs))


#combines
pheno <- pheno1 %>%
  left_join(pheno2, by="eid") %>%
  left_join(PCs %>% select(-FID), by=c("eid"="IID")) %>%
  rename(IID=eid,BMI_PGS = `26216-0.0`, BMI=`21001-0.0`) %>% drop_na()
pheno$BMI_z <- (pheno$BMI - mean(pheno$BMI)) / sd(pheno$BMI)

loc_pheno <- paste0(dir_scratch,"pheno.txt")
#fwrite(pheno, loc_pheno, sep="\t", na="NA")
pheno <- as_tibble(fread(loc_pheno))
# plots general PGS-phenotype relationship
ggplot(pheno, aes(x=BMI_z, BMI_PGS)) +
  geom_point(alpha=0.005) +
  geom_smooth(method='lm')

(cor1 <- cor.test(pheno$BMI, pheno$BMI_PGS))
cor1$estimate^2

#PCs
ggplot(pheno, aes(x=pc1, pc2)) +
  geom_point(alpha=0.005) +
  geom_smooth(method='lm')

# frNN
col_dims <- paste0("pc",1:40)
R <- 25
N <- 10000
i_keep <- sample(1:nrow(pheno),N, replace = FALSE)

pheno_N <- pheno %>%
  filter(row_number() %in% i_keep)
pheno_dim <- pheno_N %>%
  select(all_of(col_dims))
t1 <- Sys.time()
NN1 <- dbscan::frNN(pheno_dim, eps=R)
Sys.time() - t1 #26s for N=10k, R=50

pheno_N$n_neighbors <- sapply(NN1$id, length)
ggplot(pheno_dim, aes(x=pc1, y=pc2)) +
  geom_point(alpha=0.05) +
  geom_point(data=pheno_dim[NN1$id[[2]],], color="red", alpha=0.2) +
  geom_point(data=pheno_dim[NN1$id[[152]],], color="blue", alpha=0.2)

cor(pheno_N$BMI, pheno_N$BMI_PGS)
cor(pheno_N$BMI[NN1$id[[2]]], pheno_N$BMI_PGS[NN1$id[[2]]])
cor(pheno_N$BMI[NN1$id[[152]]], pheno_N$BMI_PGS[NN1$id[[152]]])

# calculates r within neighborhood
for (i in 1:nrow(pheno_N)) {
  r <- cor(pheno_N$BMI[NN1$id[[i]]], pheno_N$BMI_PGS[NN1$id[[i]]])
  pheno_N$BMI_PGS_r[i] <- r
}
for (i in 1:nrow(pheno_N)) {
  mean_r <- mean(pheno_N$BMI_PGS_r[NN1$id[[i]]])
  pheno_N$BMI_PGS_r_mean[i] <- mean_r
}
coord_cols <- paste0("pc", 1:40)
coord_means <- colMeans(select(pheno_N, all_of(coord_cols)))
pheno_N <- pheno_N %>%
  mutate(PC_dist = sqrt(rowSums((select(., all_of(coord_cols)) - coord_means)^2)))
ggplot(pheno_N[pheno_N$n_neighbors > 20,], aes(x=pc1, y=pc2)) +
  geom_point(aes(color=BMI_PGS_r_mean),alpha=0.2)

ggplot(pheno_N[pheno_N$n_neighbors > 20,],
       aes(x=PC_dist, y=BMI_PGS_r)) +
  geom_point(alpha=0.1) +
  geom_smooth(method = 'lm')

ggplot(pheno_N[pheno_N$n_neighbors > 20,],
       aes(x=PC_dist, y=BMI_PGS_r_mean^2)) +
  geom_point(alpha=0.1) +
  geom_smooth(method = 'lm')


######
IIDs_to_keep <- pheno_N$IID
loc_phenoEC <- "~/group_nuno/PXS-pipeline/scratch/pheno_EC.txt"
phenoEC <- as_tibble(fread(loc_phenoEC))
EC_i_keep <- (1:nrow(phenoEC))[phenoEC$IID %in% IIDs_to_keep]
col_expos <- colnames(phenoEC %>% select(starts_with("f",ignore.case = FALSE)))

pheno_expos <- phenoEC[EC_i_keep,c("IID",col_expos)]
#pheno_N <- pheno_N %>% left_join(phenoEC %>% select(IID, any_of(col_expos)), by="IID")

pheno_matrix <- pheno_expos %>% select(-IID)
dineof <- dineof(as.matrix(pheno_matrix))

pheno_imputed <- as_tibble(dineof$Xa)

pca <- prcomp(pheno_imputed)
pca_tibble <- as_tibble(pca$x)
colnames(pca_tibble) <- paste0("e",colnames(pca_tibble))
summary(pca)

pheno2 <- as_tibble(cbind(pheno_expos[,"IID"],pheno_imputed,pca_tibble))
pheno2 <- pheno2 %>% left_join(pheno_N, by="IID")

# plot ePC
ggplot(pheno2, aes(x=ePC1, y=ePC2)) +
  geom_point(alpha=0.1)
# compute NN again
col_dims <- paste0("ePC",1:139)
coord_means <- colMeans(select(pheno2, all_of(col_dims)))
pheno2 <- pheno2 %>%
  mutate(ePC_dist = sqrt(rowSums((select(., all_of(col_dims)) - coord_means)^2)))
R <- 9
pheno_dim <- pheno2 %>% select(all_of(col_dims))
t1 <- Sys.time()
NN2 <- dbscan::frNN(pheno_dim, eps=R)
Sys.time() - t1

pheno2$n_neighbors_ePC <- sapply(NN2$id, length)
hist(pheno2$n_neighbors_ePC)
summary(pheno2$n_neighbors_ePC)
# calculates r within neighborhood
for (i in 1:nrow(pheno2)) {
  r <- cor(pheno2$BMI[NN2$id[[i]]], pheno2$BMI_PGS[NN2$id[[i]]])
  pheno2$BMI_PGS_r_ePC[i] <- r
}
for (i in 1:nrow(pheno2)) {
  mean_r <- mean(pheno2$BMI_PGS_r[NN2$id[[i]]])
  pheno2$BMI_PGS_r_ePC_mean[i] <- mean_r
}

ggplot(pheno2, aes(x=ePC1, y=ePC2)) +
  geom_point(aes(color=BMI_PGS_r_ePC_mean),alpha=0.2)

ggplot(pheno2[pheno2$n_neighbors_ePC > 30, ],
       aes(x=ePC_dist, y=BMI_PGS_r_ePC)) +
  geom_point(alpha=0.1) +
  geom_smooth(method = 'lm')

ggplot(pheno2[which(pheno2$n_neighbors_ePC > 30,
                     pheno2$BMI_PGS_r_ePC_mean >= 0), ],
       aes(x=ePC_dist, y=BMI_PGS_r_ePC_mean^2)) +
  geom_point(alpha=0.1) +
  geom_smooth(method = 'lm') + ylim(0,0.1)
pheno2_clean  <- pheno2[which(pheno2$n_neighbors_ePC > 30,
                              pheno2$BMI_PGS_r_ePC_mean >= 0), ]
lm1 <- lm(BMI_PGS_r_ePC_mean^2 ~ ePC_dist, data=pheno2_clean)
summary(lm1)
#
ggplot(pheno2, aes(x=PC_dist, y=ePC_dist)) +
  geom_point(alpha=0.1) +
  geom_smooth(method='lm')
