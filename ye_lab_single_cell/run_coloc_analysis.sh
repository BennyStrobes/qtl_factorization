#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --nodes=1



surge_interaction_eqtl_dir="$1"
static_eqtl_dir="$2"
coloc_input_dir="$3"
coloc_dir="$4"
visualize_coloc_dir="$5"

echo $static_eqtl_dir
echo $visualize_coloc_dir

processed_gwas_studies_file=$coloc_dir"processed_gwas_studies.txt"
if false; then
python prepare_gwas_data_for_coloc_analysis.py $coloc_input_dir $processed_gwas_studies_file $coloc_dir
fi



num_components="6"
if false; then
for latent_factor_num in $(seq 1 $(($num_components))); do 
	eqtl_study_name="surge_latent_factor_"$latent_factor_num"_interaction"
	sig_eqtl_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_v3_perm_False_interaction_eqtl_results_latent_factor_"$latent_factor_num"_genome_wide_signficant_bf_fdr_0.05.txt"
	all_eqtl_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_v3_perm_False_interaction_eqtl_results_latent_factor_"$latent_factor_num"_merged.txt"
	while read gwas_study_name gwas_study_file_root; do
		echo $gwas_study_name
		echo $gwas_study_file_root
		output_root=$coloc_dir$eqtl_study_name"_"$gwas_study_name"_coloc_"
		sh run_coloc_for_study_pairing.sh $sig_eqtl_file $all_eqtl_file $gwas_study_file_root $output_root $eqtl_study_name $gwas_study_name "version1"
	done<$processed_gwas_studies_file
done
fi





eqtl_study_name="standard_eqtl"
sig_eqtl_file=$static_eqtl_dir"standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_results_genome_wide_signficant_bf_fdr_0.05.txt"
all_eqtl_file=$static_eqtl_dir"standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_results_merged.txt"
if false; then
while read gwas_study_name gwas_study_file_root; do
	echo $gwas_study_name
	output_root=$coloc_dir$eqtl_study_name"_"$gwas_study_name"_coloc_"
	sh run_coloc_for_study_pairing.sh $sig_eqtl_file $all_eqtl_file $gwas_study_file_root $output_root $eqtl_study_name $gwas_study_name "version1"
done<$processed_gwas_studies_file
fi


module load r/3.6.3
Rscript visualize_coloc_results.R $processed_gwas_studies_file $coloc_dir $visualize_coloc_dir








