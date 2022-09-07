#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=20GB





ldsc_source_code_dir="$1"
sldsc_input_data_dir="$2"
ukbb_sumstats_dir="$3"
sldsc_processed_data_dir="$4"
sldsc_results_dir="$5"
sldsc_visualization_dir="$6"
surge_interaction_eqtl_dir="$7"
sumstats_dir="$8"





################################
# Generate various versions of annotations files for SLDSC
################################
# File to be created that contains a list of all of the GWAS studies that have been processed along with the location of the processed gwas studies files
processed_gwas_studies_file=$sldsc_processed_data_dir"processed_gwas_studies.txt"
if false; then
sh prepare_gwas_data_for_sldsc_analysis.sh $ukbb_sumstats_dir $processed_gwas_studies_file $sldsc_processed_data_dir $ldsc_source_code_dir $sldsc_input_data_dir $sumstats_dir
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
if false; then
while read study_name study_file; do
	# BaselineLD eQTL
	sbatch sldsc_in_single_study.sh $study_name $study_file $ldsc_source_code_dir $sldsc_processed_data_dir $sldsc_results_dir $sldsc_input_data_dir
done <$processed_gwas_studies_file
fi



################################
# Visualize S-LDSC results
################################
if false; then
module load R/3.5.1
Rscript visualize_sldsc_results.R $processed_gwas_studies_file $sldsc_results_dir $sldsc_visualization_dir
fi



