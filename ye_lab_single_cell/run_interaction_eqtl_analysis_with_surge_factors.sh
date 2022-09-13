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
qtl_genotype_file=$eqtl_input_dir"genotype.txt"
qtl_new_covariate_file=$eqtl_input_dir"all_covariates.txt"
qtl_test_names_file=$eqtl_input_dir"variant_gene_pairs.txt"
qtl_sample_overlap_file=$eqtl_input_dir"sample_overlap.txt"



# Surge results
surge_results_stem=$eqtl_factorization_results_dir"eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_surge_results_k_10_seed_1_warm_5_rv_std_True_perm_"
surge_results_suffix="_delta_elbo_1e-2_filter_hwe_alt_init_"
surge_final_results_suffix="_delta_elbo_1e-2_filter_hwe_alt_init_"


# Extract SURGE Latent factors
surge_latent_factor_file=$output_stem"perm_False_surge_latent_factors.txt"
perm_surge_latent_factor_factor_file=$output_stem"perm_interaction_only_surge_latent_factors.txt"
if false; then
python extract_surge_latent_factor_by_pve.py $surge_results_stem $surge_latent_factor_file $perm_surge_latent_factor_factor_file $surge_results_suffix
fi

# Interaction analysis on real data
num_jobs="100"
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"perm_False_interaction_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_new_covariate_file $surge_latent_factor_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
done
fi


# Interaction analysis on permuted data
sample_permutation_file=$surge_results_stem"interaction_only"$surge_final_results_suffix"sample_permutation.txt"
num_jobs="100"
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"perm_interaction_only_interaction_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_permuted_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_new_covariate_file $perm_surge_latent_factor_factor_file $qtl_sample_overlap_file $qtl_output_root $sample_permutation_file $job_number $num_jobs
done
fi

if false; then
qtl_output_root=$output_stem"perm_False_interaction_eqtl_results_"
python2 merge_parallelized_latent_factor_interaction_eqtl_calls.py $qtl_output_root $num_jobs


qtl_output_root=$output_stem"perm_interaction_only_interaction_eqtl_results_"
python2 merge_parallelized_latent_factor_interaction_eqtl_calls.py $qtl_output_root $num_jobs
fi


module load r/3.6.3
Rscript visualize_surge_interaction_eqtls.R $output_stem"perm_"


if false; then
python2 surge_interaction_eqtl_debugger.py $output_stem"perm_"
fi


















# Interaction analysis on real data using huber sandwhich estimators
num_jobs="100"
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"perm_False_interaction_huber_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_latent_factor_interaction_huber_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_new_covariate_file $surge_latent_factor_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
done
fi


# Interaction analysis on permuted data using huber sandwhich estimators
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"perm_interaction_only_interaction_huber_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_permuted_latent_factor_interaction_huber_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_new_covariate_file $perm_surge_latent_factor_factor_file $qtl_sample_overlap_file $qtl_output_root $sample_permutation_file $job_number $num_jobs
done
fi









































##########################
# OLD
##########################


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
python merge_parallelized_latent_factor_interaction_eqtl_calls.py $output_stem"interaction_eqtl_results_v2_" $num_jobs
fi


# Interaction eqtl analysis (w standardized genotype)
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"interaction_eqtl_standardized_genotype_results_v2_"$job_number"_"$num_jobs"_"
	parallelized_expression_file=$output_stem"expression_parallelized_"$job_number".txt"
	parallelized_genotype_file=$output_stem"genotype_parallelized_"$job_number".txt"
	parallelized_test_names_file=$output_stem"variant_gene_pairs_parallelized_"$job_number".txt"
	sbatch run_latent_factor_interaction_eqtl_standardized_genotype_analysis_in_parallel.sh $parallelized_expression_file $parallelized_genotype_file $parallelized_test_names_file $qtl_new_covariate_file $surge_latent_factor_file $qtl_sample_overlap_file $qtl_output_root
done
fi

if false; then
python merge_parallelized_latent_factor_interaction_eqtl_calls_standardized_genotype.py $output_stem"interaction_eqtl_standardized_genotype_results_v2_" $num_jobs
fi









# Interaction eqtl analysis (w standardized genotype) only 3 lfs
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"interaction_eqtl_standardized_genotype_results_3_surge_contexts_"$job_number"_"$num_jobs"_"
	parallelized_expression_file=$output_stem"expression_parallelized_"$job_number".txt"
	parallelized_genotype_file=$output_stem"genotype_parallelized_"$job_number".txt"
	parallelized_test_names_file=$output_stem"variant_gene_pairs_parallelized_"$job_number".txt"
	sbatch run_latent_factor_interaction_eqtl_standardized_genotype_analysis_in_parallel.sh $parallelized_expression_file $parallelized_genotype_file $parallelized_test_names_file $qtl_new_covariate_file $surge_latent_factor_file $qtl_sample_overlap_file $qtl_output_root "3"
done
fi
if false; then
python merge_parallelized_latent_factor_interaction_eqtl_calls_standardized_genotype.py $output_stem"interaction_eqtl_standardized_genotype_results_3_surge_contexts_" $num_jobs
fi

if false; then
surge_latent_factor_file_3_components=$output_stem"surge_latent_factors_3_components_v2.txt"
cut -d$'\t' -f 1-3 $surge_latent_factor_file >$surge_latent_factor_file_3_components
fi







#########
# Perm run 1 (fixed effects and interaction)
##########
# Get pve ordered latent contexts
surge_latent_factor_file=$output_stem"perm_fi_surge_latent_factors.txt"
surge_latent_factors_file=$eqtl_factorization_results_dir"u"
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
	sbatch run_permuted_latent_factor_interaction_eqtl_analysis_in_parallel.sh $parallelized_expression_file $parallelized_genotype_file $parallelized_test_names_file $qtl_new_covariate_file $surge_latent_factor_file $qtl_sample_overlap_file $sample_permutation_file $qtl_output_root "fixed_and_interaction"
done
fi

if false; then
python merge_parallelized_latent_factor_interaction_eqtl_calls.py $output_stem"perm_fi_interaction_eqtl_results_" $num_jobs
fi




#########
# Perm run 2 (interaction only)
##########
# Get pve ordered latent contexts
surge_latent_factor_file=$output_stem"perm_i_only_surge_latent_factors.txt"
surge_latent_factors_file=$eqtl_factorization_results_dir"eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_factorization_vi_ard_heteroskedastic_results_k_10_seed_1_warm_3000_rv_std_True_perm_interaction_only_round_gen_True_temper_U_S.txt"
factor_pve_file=$eqtl_factorization_results_dir"eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_factorization_vi_ard_heteroskedastic_results_k_10_seed_1_warm_3000_rv_std_True_perm_interaction_only_round_gen_True_temper_factor_pve.txt"
if false; then
python extract_surge_latent_factor_by_pve.py $surge_latent_factors_file $factor_pve_file $surge_latent_factor_file
fi

sample_permutation_file=$eqtl_factorization_results_dir"eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_factorization_vi_ard_heteroskedastic_results_k_10_seed_1_warm_3000_rv_std_True_perm_interaction_only_round_gen_True_sample_permutation.txt"
# Interaction eqtl analysis 2
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"perm_i_only_interaction_eqtl_results_"$job_number"_"$num_jobs"_"
	parallelized_expression_file=$output_stem"expression_parallelized_"$job_number".txt"
	parallelized_genotype_file=$output_stem"genotype_parallelized_"$job_number".txt"
	parallelized_test_names_file=$output_stem"variant_gene_pairs_parallelized_"$job_number".txt"
	sbatch run_permuted_latent_factor_interaction_eqtl_analysis_in_parallel.sh $parallelized_expression_file $parallelized_genotype_file $parallelized_test_names_file $qtl_new_covariate_file $surge_latent_factor_file $qtl_sample_overlap_file $sample_permutation_file $qtl_output_root "interaction_only"
done
fi




if false; then
module load R/3.5.1
Rscript visualize_surge_interaction_eqtls.R $output_stem
fi





# Standard eqtl analysis with standardized genotype
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"standard_eqtl_standardized_genotype_results_"$job_number"_"$num_jobs"_"
	parallelized_expression_file=$output_stem"expression_parallelized_"$job_number".txt"
	parallelized_genotype_file=$output_stem"genotype_parallelized_"$job_number".txt"
	parallelized_test_names_file=$output_stem"variant_gene_pairs_parallelized_"$job_number".txt"
	sbatch run_standard_eqtl_standardized_genotype_analysis_in_parallel.sh $parallelized_expression_file $parallelized_genotype_file $parallelized_test_names_file $qtl_new_covariate_file $qtl_sample_overlap_file $qtl_output_root
done
fi

if false; then
python merge_parallelized_latent_factor_interaction_eqtl_calls_standardized_genotype.py $output_stem"standard_eqtl_standardized_genotype_results_" $num_jobs
fi








