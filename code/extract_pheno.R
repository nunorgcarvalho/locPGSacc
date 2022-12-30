library(tidyverse)
library(data.table)

# loads paths, directories, and definitions
source("paths_and_definitions.R")

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
ancestry_table <- tibble(FID = PCs$FID, IID = PCs$IID, ancestry = group)
loc_out <- "../scratch/ancestry_table.txt"
write.table(ancestry_table, loc_out, sep="\t", quote=FALSE, row.names = FALSE)
