#!/bin/bash

dir_repo="/n/groups/patel/nuno/locPGSacc/"
loc_cols=${dir_repo}'input_data/ukb_pheno_cols.txt'
loc_pheno1_full="/n/groups/patel/uk_biobank/project_22881_669542/ukb669542.csv"
#loc_pheno2_full="/n/groups/patel/uk_biobank/main_data_34521/ukb34521.csv"
loc_pheno2_full=${dir_repo}'scratch/ukb34521.csv'
loc_pheno1_filtered=${dir_repo}'scratch/pheno1_filtered.csv'
loc_pheno2_filtered=${dir_repo}'scratch/pheno2_filtered.csv'

# Read columns from $loc_cols into an array
mapfile -t cols < "$loc_cols"

# Convert the array to a comma-separated string
cols_csv=$(IFS=,; echo "${cols[*]}")

# Generate the awk command
awk_command=$(printf 'BEGIN {
FS="\",\""
OFS=","
} 
NR==1 {
    for (i=1; i<=NF; i++) {
        f[$i] = i
    }
    split("%s", a, ",")
}
{
    for (i in a) {
        if (a[i] in f) {
            printf "%s%s" $(f[a[i]]), OFS
        }
    }
    print ""
}' "$cols_csv")

# Apply awk command to filter columns of $loc_pheno2_full and save to $loc_pheno2_filtered
awk "$awk_command" "$loc_pheno2_full" > "$loc_pheno2_filtered"
