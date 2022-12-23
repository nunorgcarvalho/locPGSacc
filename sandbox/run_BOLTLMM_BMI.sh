#!/bin/sh
#SBATCH -c 20
#SBATCH -t 2-00:00
#SBATCH -p medium
#SBATCH --mem=125G
#SBATCH -o BMI_BOLTLMM.out
#SBATCH -e BMI_BOLTLMM.err

dir_scratch="~/scratch3/locPGSacc/"

~/bolt \
--numThreads 20 \
--bed /n/groups/patel/uk_biobank/main_data_9512/ukb_cal_chr{1:22}_v2.bed \
--bim /n/groups/patel/uk_biobank/main_data_9512/ukb_snp_chr{1:22}_v2.bim \
--fam /n/groups/patel/uk_biobank/main_data_9512/ukb_bolt_lmm.fam \
--LDscoresFile /n/groups/patel/bin/BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz \
--remove ~/jobs/PXS_pipeline/input_data/bolt.in_plink_but_not_imputed.FID_IID.978.txt \
--phenoFile ${dir_scratch}pheno_BMI_train.txt \
--phenoCol BMI \
--covarFile ${dir_scratch}pheno_BMI_train.txt \
--covarCol sex \
--covarCol assessment_center \
--qCovarCol age \
--qCovarCol pc{1:40} \
--covarMaxLevels 25 \
--lmm \
--verboseStats \
--statsFile ${dir_scratch}LMM_BMI.txt \
--bgenFile /n/no_backup2/patel/ukb_imp_chr{1:22}_v3.bgen \
--bgenMinMAF 1e-3 \
--bgenMinINFO 0.3 \
--sampleFile /n/no_backup2/patel/ukb22881_imp_chr1_v3_s487324.sample \
--statsFileBgenSnps ${dir_scratch}LMM_BMI_bgen.txt
