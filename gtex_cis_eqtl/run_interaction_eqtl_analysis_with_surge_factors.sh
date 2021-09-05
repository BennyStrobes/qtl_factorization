#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1



eqtl_input_dir="$1"
surge_latent_factors_file="$2"
factor_pve_file="$3"
latent_factor_num="$4"
output_stem="$5"


# Input data
qtl_expression_file=$eqtl_input_dir"cross_tissue_eqtl_expression_input.txt"
qtl_genotype_file=$eqtl_input_dir"cross_tissue_eqtl_genotype_input.txt"
qtl_covariate_file=$eqtl_input_dir"cross_tissue_eqtl_covariate_input.txt"
qtl_test_names_file=$eqtl_input_dir"all_tests.txt"
qtl_sample_overlap_file=$eqtl_input_dir"individual_id.txt"


surge_latent_factor_file=$output_stem"surge_latent_factor.txt"
if false; then
python extract_surge_latent_factor_by_pve.py $surge_latent_factors_file $factor_pve_file $latent_factor_num $surge_latent_factor_file
fi

num_jobs="1"
job_number="0"
qtl_output_root=$output_stem"interaction_eqtl_results_"$job_number"_"$num_jobs"_"
if false; then
sbatch run_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $surge_latent_factor_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
fi


qtl_output_root=$output_stem"interaction_eqtl_results_"
python merge_parallelized_interaction_eqtl_calls.py $qtl_output_root $num_jobs