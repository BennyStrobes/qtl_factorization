#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --mem=20GB



coloc_input_dir="$1"
processed_gwas_studies_file="$2"
sldsc_processed_data_dir="$3"
ldsc_source_code_dir="$4"
sldsc_input_data_dir="$5"
sumstats_dir="$6"

source ~/.bash_profile
conda activate ldsc

study_name="ukbb_blood_monocyte_count"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" > $processed_gwas_studies_file
python remove_very_low_pvalues_from_ukbb_summary_statistics_file.py ${coloc_input_dir}blood_MONOCYTE_COUNT.sumstats.gz ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz"
python ${ldsc_source_code_dir}munge_sumstats.py --sumstats ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz" --merge-alleles ${sldsc_input_data_dir}w_hm3.snplist --out ${sldsc_processed_data_dir}${study_name} --a1-inc


study_name="ukbb_blood_lymphocyte_count"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
python remove_very_low_pvalues_from_ukbb_summary_statistics_file.py ${coloc_input_dir}blood_LYMPHOCYTE_COUNT.sumstats.gz ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz"
python ${ldsc_source_code_dir}munge_sumstats.py --sumstats ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz" --merge-alleles ${sldsc_input_data_dir}w_hm3.snplist --out ${sldsc_processed_data_dir}${study_name} --a1-inc



study_name="ukbb_eczema"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
python remove_very_low_pvalues_from_ukbb_summary_statistics_file.py ${coloc_input_dir}disease_ALLERGY_ECZEMA_DIAGNOSED.sumstats.gz ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz"
python ${ldsc_source_code_dir}munge_sumstats.py --sumstats ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz" --merge-alleles ${sldsc_input_data_dir}w_hm3.snplist --out ${sldsc_processed_data_dir}${study_name} --a1-inc


study_name="ukbb_bmi"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
python remove_very_low_pvalues_from_ukbb_summary_statistics_file.py ${coloc_input_dir}body_BMIz.sumstats.gz ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz"
python ${ldsc_source_code_dir}munge_sumstats.py --sumstats ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz" --merge-alleles ${sldsc_input_data_dir}w_hm3.snplist --out ${sldsc_processed_data_dir}${study_name} --a1-inc



study_name="ukbb_blood_eosinophil_count"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
python remove_very_low_pvalues_from_ukbb_summary_statistics_file.py ${coloc_input_dir}blood_EOSINOPHIL_COUNT.sumstats.gz ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz"
python ${ldsc_source_code_dir}munge_sumstats.py --sumstats ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz" --merge-alleles ${sldsc_input_data_dir}w_hm3.snplist --out ${sldsc_processed_data_dir}${study_name} --a1-inc



study_name="ukbb_blood_high_light_scatter_reticulotye_count"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
python remove_very_low_pvalues_from_ukbb_summary_statistics_file.py ${coloc_input_dir}blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT.sumstats.gz ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz"
python ${ldsc_source_code_dir}munge_sumstats.py --sumstats ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz" --merge-alleles ${sldsc_input_data_dir}w_hm3.snplist --out ${sldsc_processed_data_dir}${study_name} --a1-inc


study_name="ukbb_blood_mean_corpuscular_hemoglobin"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
python remove_very_low_pvalues_from_ukbb_summary_statistics_file.py ${coloc_input_dir}blood_MEAN_CORPUSCULAR_HEMOGLOBIN.sumstats.gz ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz"
python ${ldsc_source_code_dir}munge_sumstats.py --sumstats ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz" --merge-alleles ${sldsc_input_data_dir}w_hm3.snplist --out ${sldsc_processed_data_dir}${study_name} --a1-inc



study_name="ukbb_blood_platelet_vol"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
python remove_very_low_pvalues_from_ukbb_summary_statistics_file.py ${coloc_input_dir}blood_MEAN_PLATELET_VOL.sumstats.gz ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz"
python ${ldsc_source_code_dir}munge_sumstats.py --sumstats ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz" --merge-alleles ${sldsc_input_data_dir}w_hm3.snplist --out ${sldsc_processed_data_dir}${study_name} --a1-inc



study_name="ukbb_blood_platelet_count"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
python remove_very_low_pvalues_from_ukbb_summary_statistics_file.py ${coloc_input_dir}blood_PLATELET_COUNT.sumstats.gz ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz"
python ${ldsc_source_code_dir}munge_sumstats.py --sumstats ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz" --merge-alleles ${sldsc_input_data_dir}w_hm3.snplist --out ${sldsc_processed_data_dir}${study_name} --a1-inc


study_name="ukbb_blood_red_count"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
python remove_very_low_pvalues_from_ukbb_summary_statistics_file.py ${coloc_input_dir}blood_RED_COUNT.sumstats.gz ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz"
python ${ldsc_source_code_dir}munge_sumstats.py --sumstats ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz" --merge-alleles ${sldsc_input_data_dir}w_hm3.snplist --out ${sldsc_processed_data_dir}${study_name} --a1-inc


study_name="ukbb_blood_white_count"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
python remove_very_low_pvalues_from_ukbb_summary_statistics_file.py ${coloc_input_dir}blood_WHITE_COUNT.sumstats.gz ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz"
python ${ldsc_source_code_dir}munge_sumstats.py --sumstats ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz" --merge-alleles ${sldsc_input_data_dir}w_hm3.snplist --out ${sldsc_processed_data_dir}${study_name} --a1-inc



study_name="ukbb_height"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
python remove_very_low_pvalues_from_ukbb_summary_statistics_file.py ${coloc_input_dir}body_HEIGHTz.sumstats.gz ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz"
python ${ldsc_source_code_dir}munge_sumstats.py --sumstats ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz" --merge-alleles ${sldsc_input_data_dir}w_hm3.snplist --out ${sldsc_processed_data_dir}${study_name} --a1-inc


study_name="ukbb_T2D"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
python remove_very_low_pvalues_from_ukbb_summary_statistics_file.py ${coloc_input_dir}disease_T2D.sumstats.gz ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz"
python ${ldsc_source_code_dir}munge_sumstats.py --sumstats ${sldsc_processed_data_dir}${study_name}"_low_pvalue_removed.sumstats.gz" --merge-alleles ${sldsc_input_data_dir}w_hm3.snplist --out ${sldsc_processed_data_dir}${study_name} --a1-inc



study_name="Celiac"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
cp ${sumstats_dir}PASS_Celiac.sumstats ${sldsc_processed_data_dir}${study_name}".sumstats"
gzip ${sldsc_processed_data_dir}${study_name}".sumstats"

study_name="Crohns"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
cp ${sumstats_dir}PASS_Crohns_Disease.sumstats ${sldsc_processed_data_dir}${study_name}".sumstats"
gzip ${sldsc_processed_data_dir}${study_name}".sumstats"


study_name="Ulcerative_Colitis"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
cp ${sumstats_dir}PASS_Ulcerative_Colitis.sumstats ${sldsc_processed_data_dir}${study_name}".sumstats"
gzip ${sldsc_processed_data_dir}${study_name}".sumstats"


study_name="Rheumatoid_Arthritis"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
cp ${sumstats_dir}PASS_Rheumatoid_Arthritis.sumstats ${sldsc_processed_data_dir}${study_name}".sumstats"
gzip ${sldsc_processed_data_dir}${study_name}".sumstats"


study_name="Lupus"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
cp ${sumstats_dir}PASS_Lupus.sumstats ${sldsc_processed_data_dir}${study_name}".sumstats"
gzip ${sldsc_processed_data_dir}${study_name}".sumstats"



study_name="IBD"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
cp ${sumstats_dir}PASS_IBD.sumstats ${sldsc_processed_data_dir}${study_name}".sumstats"
gzip ${sldsc_processed_data_dir}${study_name}".sumstats"


study_name="Multiple_sclerosis"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
cp ${sumstats_dir}PASS_Multiple_sclerosis.sumstats ${sldsc_processed_data_dir}${study_name}".sumstats"
gzip ${sldsc_processed_data_dir}${study_name}".sumstats"



study_name="PBC"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
cp ${sumstats_dir}PASS_Primary_biliary_cirrhosis.sumstats ${sldsc_processed_data_dir}${study_name}".sumstats"
gzip ${sldsc_processed_data_dir}${study_name}".sumstats"



study_name="CAD"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
cp ${sumstats_dir}PASS_Coronary_Artery_Disease.sumstats ${sldsc_processed_data_dir}${study_name}".sumstats"
gzip ${sldsc_processed_data_dir}${study_name}".sumstats"


study_name="Bipolar"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
cp ${sumstats_dir}PASS_Bipolar_Disorder.sumstats ${sldsc_processed_data_dir}${study_name}".sumstats"
gzip ${sldsc_processed_data_dir}${study_name}".sumstats"


study_name="Alzheimer"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
cp ${sumstats_dir}PASS_Alzheimer.sumstats ${sldsc_processed_data_dir}${study_name}".sumstats"
gzip ${sldsc_processed_data_dir}${study_name}".sumstats"


study_name="Schizophrenia"
echo -e "$study_name\t"${sldsc_processed_data_dir}${study_name}".sumstats.gz" >> $processed_gwas_studies_file
cp ${sumstats_dir}PASS_Schizophrenia.sumstats ${sldsc_processed_data_dir}${study_name}".sumstats"
gzip ${sldsc_processed_data_dir}${study_name}".sumstats"




