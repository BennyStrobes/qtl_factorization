#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


surge_interaction_eqtl_dir="$1"
processed_genotype_dir="$2"
susie_input_data_dir="$3"
susie_results_dir="$4"
susie_visualization_dir="$5"
sample_names_file="$6"



# Input data
latent_factor_num="1"
eqtl_study_name="surge_latent_factor_"$latent_factor_num"_interaction"
sig_eqtl_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_interaction_eqtl_results_latent_factor_"$latent_factor_num"_genome_wide_signficant_bf_fdr_0.05.txt"
all_eqtl_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_interaction_eqtl_results_latent_factor_"$latent_factor_num"_merged.txt"
surge_latent_factor_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_surge_latent_factors.txt"
# Output files
output_processed_data_root=$susie_input_data_dir$eqtl_study_name"_susie_finemapping_"
output_results_root=$susie_results_dir$eqtl_study_name"_susie_finemapping_"
output_visualization_root=$susie_visualization_dir$eqtl_study_name"_susie_finemapping_"

sh run_susie_finemapping_in_single_study.sh $sig_eqtl_file $all_eqtl_file $eqtl_study_name $processed_genotype_dir $output_processed_data_root $output_results_root $output_visualization_root $surge_latent_factor_file $latent_factor_num $sample_names_file