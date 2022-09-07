#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --mem=20GB




source ~/.bash_profile
conda activate ldsc



chrom_num="$1"
sldsc_input_data_dir="$2"
sldsc_processed_data_dir="$3"
ldsc_source_code_dir="$4"

cp ${sldsc_input_data_dir}"1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD."$chrom_num".annot.gz" ${sldsc_processed_data_dir}
python ${ldsc_source_code_dir}ldsc.py --l2 --bfile ${sldsc_input_data_dir}"1000G_EUR_Phase3_plink/1000G.EUR.QC."${chrom_num} --ld-wind-cm 1 --annot ${sldsc_input_data_dir}"1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD."$chrom_num".annot.gz" --out ${sldsc_processed_data_dir}"baselineLD."$chrom_num --print-snps ${sldsc_input_data_dir}"hapmap3_snps/hm."${chrom_num}".snp"