#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=10GB


study_name="$1"
study_file="$2"
ldsc_source_code_dir="$3"
sldsc_processed_data_dir="$4"
sldsc_results_dir="$5"
sldsc_input_data_dir="$6"

module load python/2.7-anaconda


echo "SLDSC in "$study_name
if false; then
python ${ldsc_source_code_dir}ldsc.py --h2 ${study_file} --ref-ld-chr ${sldsc_processed_data_dir}"baselineLD.,"${sldsc_processed_data_dir}"surge_eqtls_1e-05." --w-ld-chr ${sldsc_input_data_dir}"1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${sldsc_input_data_dir}"1000G_Phase3_frq/1000G.EUR.QC." --out ${sldsc_results_dir}${study_name}"_surge_eqtls_1e-05_and_baselineLD"
	
python ${ldsc_source_code_dir}ldsc.py --h2 ${study_file} --ref-ld-chr ${sldsc_processed_data_dir}"baselineLD.,"${sldsc_processed_data_dir}"joint_surge_eqtls_1e-05." --w-ld-chr ${sldsc_input_data_dir}"1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${sldsc_input_data_dir}"1000G_Phase3_frq/1000G.EUR.QC." --out ${sldsc_results_dir}${study_name}"_joint_surge_eqtls_1e-05_and_baselineLD"
fi

# Surge eqtls in each component seperately
for component_num in $(seq 1 $((10))); do 
	echo $component_num
	python ${ldsc_source_code_dir}ldsc.py --h2 ${study_file} --ref-ld-chr ${sldsc_processed_data_dir}"baselineLD.,"${sldsc_processed_data_dir}"surge_"${component_num}"_eqtls_1e-05." --w-ld-chr ${sldsc_input_data_dir}"1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${sldsc_input_data_dir}"1000G_Phase3_frq/1000G.EUR.QC." --out ${sldsc_results_dir}${study_name}"_surge_"${component_num}"_eqtls_1e-05_only_and_baselineLD"
done

# Standard eqtls
python ${ldsc_source_code_dir}ldsc.py --h2 ${study_file} --ref-ld-chr ${sldsc_processed_data_dir}"baselineLD.,"${sldsc_processed_data_dir}"standard_eqtls_1e-05." --w-ld-chr ${sldsc_input_data_dir}"1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${sldsc_input_data_dir}"1000G_Phase3_frq/1000G.EUR.QC." --out ${sldsc_results_dir}${study_name}"_standard_eqtls_1e-05_only_and_baselineLD"
