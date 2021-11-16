#!/bin/bash -l

#SBATCH
#SBATCH --time=30:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=20GB



eqtl_input_dir="$1"
eqtl_factorization_results_dir="$2"
output_stem="$3"


# Input data
qtl_expression_file=$eqtl_input_dir"expression.txt"
qtl_genotype_file=$eqtl_input_dir"genotype2.txt"
qtl_covariate_file=$eqtl_input_dir"covariates.txt"
qtl_lf_file=$eqtl_input_dir"latent_factors.txt"
qtl_test_names_file=$eqtl_input_dir"variant_gene_pairs.txt"
qtl_sample_overlap_file=$eqtl_input_dir"sample_overlap.txt"


# Generate organized covariate files
qtl_new_covariate_file=$output_stem"merged_covariates.txt"
if false; then
python merge_covariates.py $qtl_covariate_file $qtl_lf_file $qtl_new_covariate_file
fi


# Split up expression and genotype file for parallelelization purposes
num_jobs="200"
total_lines=`wc -l $qtl_test_names_file`
if false; then 
python split_file_into_X_chunks.py $qtl_genotype_file $total_lines $num_jobs $output_stem"genotype_parallelized_" "False"
python split_file_into_X_chunks.py $qtl_expression_file $total_lines $num_jobs $output_stem"expression_parallelized_" "False"
python split_file_into_X_chunks.py $qtl_test_names_file $total_lines $num_jobs $output_stem"variant_gene_pairs_parallelized_" "True"
fi


# Get pve ordered latent contexts
surge_latent_factor_file=$output_stem"surge_latent_factors_v2.txt"
surge_latent_factors_file=$eqtl_factorization_results_dir"eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_factorization_vi_ard_heteroskedastic_results_k_10_seed_1_warm_3000_rv_std_True_perm_False_round_gen_True_temper_U_S.txt"
factor_pve_file=$eqtl_factorization_results_dir"eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_factorization_vi_ard_heteroskedastic_results_k_10_seed_1_warm_3000_rv_std_True_perm_False_round_gen_True_temper_factor_pve.txt"
if false; then
python extract_surge_latent_factor_by_pve.py $surge_latent_factors_file $factor_pve_file $surge_latent_factor_file
fi



# Interaction eqtl analysis
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"interaction_eqtl_results_v2_"$job_number"_"$num_jobs"_"
	parallelized_expression_file=$output_stem"expression_parallelized_"$job_number".txt"
	parallelized_genotype_file=$output_stem"genotype_parallelized_"$job_number".txt"
	parallelized_test_names_file=$output_stem"variant_gene_pairs_parallelized_"$job_number".txt"
	sbatch run_latent_factor_interaction_eqtl_analysis_in_parallel.sh $parallelized_expression_file $parallelized_genotype_file $parallelized_test_names_file $qtl_new_covariate_file $surge_latent_factor_file $qtl_sample_overlap_file $qtl_output_root
done
fi

if false; then
python merge_parallelized_latent_factor_interaction_eqtl_calls.py $output_stem"interaction_eqtl_results_" $num_jobs
fi







#########
# Perm run
##########
# Get pve ordered latent contexts
surge_latent_factor_file=$output_stem"perm_fi_surge_latent_factors.txt"
surge_latent_factors_file=$eqtl_factorization_results_dir"eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_factorization_vi_ard_heteroskedastic_results_k_10_seed_1_warm_3000_rv_std_True_perm_fixed_and_interaction_round_gen_True_temper_U_S.txt"
factor_pve_file=$eqtl_factorization_results_dir"eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_factorization_vi_ard_heteroskedastic_results_k_10_seed_1_warm_3000_rv_std_True_perm_fixed_and_interaction_round_gen_True_temper_factor_pve.txt"
if false; then
python extract_surge_latent_factor_by_pve.py $surge_latent_factors_file $factor_pve_file $surge_latent_factor_file
fi


sample_permutation_file=$eqtl_factorization_results_dir"eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_factorization_vi_ard_heteroskedastic_results_k_10_seed_1_warm_3000_rv_std_True_perm_fixed_and_interaction_round_gen_True_t_sample_permutation.txt"
# Interaction eqtl analysis 1
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"perm_fi_interaction_eqtl_results_"$job_number"_"$num_jobs"_"
	parallelized_expression_file=$output_stem"expression_parallelized_"$job_number".txt"
	parallelized_genotype_file=$output_stem"genotype_parallelized_"$job_number".txt"
	parallelized_test_names_file=$output_stem"variant_gene_pairs_parallelized_"$job_number".txt"
	sbatch run_permuted_latent_factor_interaction_eqtl_analysis_in_parallel.sh $parallelized_expression_file $parallelized_genotype_file $parallelized_test_names_file $qtl_new_covariate_file $surge_latent_factor_file $qtl_sample_overlap_file $sample_permutation_file $qtl_output_root
done
fi

# Get pve ordered latent contexts
surge_latent_factor_file=$output_stem"perm_i_only_surge_latent_factors.txt"
surge_latent_factors_file=$eqtl_factorization_results_dir"eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_factorization_vi_ard_heteroskedastic_results_k_10_seed_1_warm_3000_rv_std_True_perm_interaction_only_round_gen_True_temper_U_S.txt"
factor_pve_file=$eqtl_factorization_results_dir"eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_factorization_vi_ard_heteroskedastic_results_k_10_seed_1_warm_3000_rv_std_True_perm_interaction_only_round_gen_True_temper_factor_pve.txt"
if false; then
python extract_surge_latent_factor_by_pve.py $surge_latent_factors_file $factor_pve_file $surge_latent_factor_file
fi


if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"perm_i_only_interaction_eqtl_results_"$job_number"_"$num_jobs"_"
	parallelized_expression_file=$output_stem"expression_parallelized_"$job_number".txt"
	parallelized_genotype_file=$output_stem"genotype_parallelized_"$job_number".txt"
	parallelized_test_names_file=$output_stem"variant_gene_pairs_parallelized_"$job_number".txt"
	sbatch run_permuted_interaction_only_latent_factor_interaction_eqtl_analysis_in_parallel.sh $parallelized_expression_file $parallelized_genotype_file $parallelized_test_names_file $qtl_new_covariate_file $surge_latent_factor_file $qtl_sample_overlap_file $sample_permutation_file $qtl_output_root
done
fi



if false; then
python merge_parallelized_latent_factor_interaction_eqtl_calls.py $output_stem"perm_fi_interaction_eqtl_results_" $num_jobs
fi














# Standard eqtl analysis
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"standard_eqtl_results_"$job_number"_"$num_jobs"_"
	parallelized_expression_file=$output_stem"expression_parallelized_"$job_number".txt"
	parallelized_genotype_file=$output_stem"genotype_parallelized_"$job_number".txt"
	parallelized_test_names_file=$output_stem"variant_gene_pairs_parallelized_"$job_number".txt"
	sbatch run_standard_eqtl_analysis_in_parallel.sh $parallelized_expression_file $parallelized_genotype_file $parallelized_test_names_file $qtl_new_covariate_file $qtl_sample_overlap_file $qtl_output_root
done
fi
if false; then
python merge_parallelized_standard_eqtl_calls.py $output_stem"standard_eqtl_results_" $num_jobs
fi

if false; then
module load R/3.5.1
Rscript temp_viz.R $output_stem
fi

