#!/bin/bash -l

#SBATCH
#SBATCH --time=30:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=20GB





ldsc_source_code_dir="$1"
sldsc_input_data_dir="$2"
coloc_input_dir="$3"
sldsc_processed_data_dir="$4"
sldsc_results_dir="$5"
sldsc_visualization_dir="$6"
surge_interaction_eqtl_dir="$7"
sumstats_dir="$8"


module load python/2.7-anaconda




################################
# Generate various versions of annotations files for SLDSC
################################
# File to be created that contains a list of all of the GWAS studies that have been processed along with the location of the processed gwas studies files
processed_gwas_studies_file=$sldsc_processed_data_dir"processed_gwas_studies.txt"
if false; then
sh prepare_gwas_data_for_sldsc_analysis.sh $coloc_input_dir $processed_gwas_studies_file $sldsc_processed_data_dir $ldsc_source_code_dir $sldsc_input_data_dir $sumstats_dir
fi




################################
# Generate various versions of annotations files for SLDSC
################################
if false; then
sh generate_surge_annotation_file.sh $sldsc_input_data_dir $surge_interaction_eqtl_dir $sldsc_processed_data_dir $ldsc_source_code_dir
fi



################################
# Run S-LDSC
################################

# Loop through GWAS studies
# for each study from s-ldsc
while read study_name study_file; do
	echo $study_name

	# Using egene annotations for each surge study independently while controlling for base line anno
	python ${ldsc_source_code_dir}ldsc.py --h2 ${study_file} --ref-ld-chr ${sldsc_processed_data_dir}"baseline.,"${sldsc_processed_data_dir}"surge_egenes." --w-ld-chr ${sldsc_input_data_dir}"1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." --overlap-annot --frqfile-chr ${sldsc_input_data_dir}"1000G_Phase3_frq/1000G.EUR.QC." --out ${sldsc_results_dir}${study_name}"_surge_egenes_05_and_baseline"

	# Using egene annotations union surge studies while controlling for base line anno
	python ${ldsc_source_code_dir}ldsc.py --h2 ${study_file} --ref-ld-chr ${sldsc_processed_data_dir}"baseline.,"${sldsc_processed_data_dir}"joint_surge_egenes_05." --w-ld-chr ${sldsc_input_data_dir}"1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." --overlap-annot --frqfile-chr ${sldsc_input_data_dir}"1000G_Phase3_frq/1000G.EUR.QC." --out ${sldsc_results_dir}${study_name}"_joint_surge_egenes_05_and_baseline"
done <$processed_gwas_studies_file










