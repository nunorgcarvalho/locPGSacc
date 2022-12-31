# sets directories and paths (do NOT include last forward slash / )
dir_repo="/home/nur479/scratch3/locPGSacc"

mkdir -p ${dir_repo}/LMM/

traits=(height BMI T2D)
for trait in ${traits[@]}
do

cd ${dir_repo}/code

########################################
## Creates BOLT-LMM script and submits##
########################################
echo '#!/bin/sh
#SBATCH -c 20
#SBATCH -t 4-23:59
#SBATCH -p medium
#SBATCH --mem=125G
#SBATCH -o BOLTLMM_'${trait}'.out
#SBATCH -e BOLTLMM_'${trait}'.err
~/bolt \
--numThreads 20 \
--bed /n/groups/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/groups/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/groups/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove '${dir_repo}'/input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile '${dir_repo}'/scratch/pheno.txt \
--phenoCol '${trait}' \
--covarFile '${dir_scratch}'/scratch/pheno.txt \
--covarCol sex \
--covarCol assessment_center \
--qCovarCol age \
--qCovarCol pc{1:40} \
--covarMaxLevels 25 \
--lmm \
--verboseStats \
--statsFile '${dir_repo}'/scratch/LMM/LMM_'${trait}'.txt \
--bgenFile /n/no_backup2/patel/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF 1e-3 \
--bgenMinINFO 0.3 \
--sampleFile /n/no_backup2/patel/ukb22881_imp_chr1_v3_s487324.sample \
--statsFileBgenSnps '${dir_repo}'/scratch/LMM/LMM_'${trait}'_bgen.txt
' > ${dir_repo}/code/run_BOLTLMM_${trait}.sh
sbatch ${dir_repo}/code/run_BOLTLMM_${trait}.sh
echo 'Submitted BOLT-LMM for '${trait}

done