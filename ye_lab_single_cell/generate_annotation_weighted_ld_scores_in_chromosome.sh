#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=20GB


ldsc_source_code_dir="$1"
sldsc_input_data_dir="$2"
sldsc_processed_data_dir="$3"
chrom_num="$4"

echo "Anno for chrom "$chrom_num
module load python/2.7-anaconda


# Joint surge eqtls
if false; then
python ${ldsc_source_code_dir}ldsc.py --l2 --bfile ${sldsc_input_data_dir}"1000G_EUR_Phase3_plink/1000G.EUR.QC."${chrom_num} --ld-wind-cm 1 --annot ${sldsc_processed_data_dir}"joint_surge_eqtls_1e-05."${chrom_num}".annot" --out ${sldsc_processed_data_dir}"joint_surge_eqtls_1e-05."${chrom_num} --print-snps ${sldsc_input_data_dir}"hapmap3_snps/hm."${chrom_num}".snp"
fi

# surge eqtls
if false; then
python ${ldsc_source_code_dir}ldsc.py --l2 --bfile ${sldsc_input_data_dir}"1000G_EUR_Phase3_plink/1000G.EUR.QC."${chrom_num} --ld-wind-cm 1 --annot ${sldsc_processed_data_dir}"surge_eqtls_1e-05."${chrom_num}".annot" --out ${sldsc_processed_data_dir}"surge_eqtls_1e-05."${chrom_num} --print-snps ${sldsc_input_data_dir}"hapmap3_snps/hm."${chrom_num}".snp"
fi

# Single component surge eqtls
for component_num in $(seq 1 $((10))); do 
	python ${ldsc_source_code_dir}ldsc.py --l2 --bfile ${sldsc_input_data_dir}"1000G_EUR_Phase3_plink/1000G.EUR.QC."${chrom_num} --ld-wind-cm 1 --annot ${sldsc_processed_data_dir}"surge_"${component_num}"_eqtls_1e-05."${chrom_num}".annot" --out ${sldsc_processed_data_dir}"surge_"${component_num}"_eqtls_1e-05."${chrom_num} --print-snps ${sldsc_input_data_dir}"hapmap3_snps/hm."${chrom_num}".snp"
done


# Standard eQTLs
python ${ldsc_source_code_dir}ldsc.py --l2 --bfile ${sldsc_input_data_dir}"1000G_EUR_Phase3_plink/1000G.EUR.QC."${chrom_num} --ld-wind-cm 1 --annot ${sldsc_processed_data_dir}"standard_eqtls_1e-05."${chrom_num}".annot" --out ${sldsc_processed_data_dir}"standard_eqtls_1e-05."${chrom_num} --print-snps ${sldsc_input_data_dir}"hapmap3_snps/hm."${chrom_num}".snp"



# Baseline
if false; then
cp ${sldsc_input_data_dir}"baseline_v1.2/baseline."$chrom_num".annot.gz" ${sldsc_processed_data_dir}
python ${ldsc_source_code_dir}ldsc.py --l2 --bfile ${sldsc_input_data_dir}"1000G_EUR_Phase3_plink/1000G.EUR.QC."${chrom_num} --ld-wind-cm 1 --annot ${sldsc_processed_data_dir}"baseline."$chrom_num".annot.gz" --out ${sldsc_processed_data_dir}"baseline."$chrom_num --print-snps ${sldsc_input_data_dir}"hapmap3_snps/hm."${chrom_num}".snp"
fi


# BaselineLD
if false; then
cp ${sldsc_input_data_dir}"1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD."$chrom_num".annot.gz" ${sldsc_processed_data_dir}
python ${ldsc_source_code_dir}ldsc.py --l2 --bfile ${sldsc_input_data_dir}"1000G_EUR_Phase3_plink/1000G.EUR.QC."${chrom_num} --ld-wind-cm 1 --annot ${sldsc_processed_data_dir}"baselineLD."$chrom_num".annot.gz" --out ${sldsc_processed_data_dir}"baselineLD."$chrom_num --print-snps ${sldsc_input_data_dir}"hapmap3_snps/hm."${chrom_num}".snp"
fi