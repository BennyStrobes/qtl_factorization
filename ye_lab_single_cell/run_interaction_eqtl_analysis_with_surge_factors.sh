#!/bin/bash -l

#SBATCH
#SBATCH --time=30:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1



eqtl_input_dir="$1"
surge_latent_factors_file="$2"
factor_pve_file="$3"
latent_factor_num="$4"
output_stem="$5"


# Input data
qtl_expression_file=$eqtl_input_dir"expression.txt"
qtl_genotype_file=$eqtl_input_dir"genotype2.txt"
qtl_covariate_file=$eqtl_input_dir"covariates.txt"
qtl_lf_file=$eqtl_input_dir"latent_factors.txt"
qtl_test_names_file=$eqtl_input_dir"variant_gene_pairs.txt"
qtl_sample_overlap_file=$eqtl_input_dir"sample_overlap.txt"

qtl_new_covariate_file=$output_stem"merged_covariates.txt"
if false; then
python merge_covariates.py $qtl_covariate_file $qtl_lf_file $qtl_new_covariate_file
fi


surge_latent_factor_file=$output_stem"surge_latent_factor.txt"
if false; then
python extract_surge_latent_factor_by_pve.py $surge_latent_factors_file $factor_pve_file $latent_factor_num $surge_latent_factor_file
fi

num_jobs="200"

if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"interaction_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_new_covariate_file $surge_latent_factor_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
done
fi

python merge_parallelized_latent_factor_interaction_eqtl_calls.py $output_stem"interaction_eqtl_results_" $num_jobs


