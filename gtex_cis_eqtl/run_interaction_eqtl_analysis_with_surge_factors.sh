#!/bin/bash -l

#SBATCH
#SBATCH --time=4:00:00
#SBATCH --nodes=1




eqtl_input_dir="$1"
surge_results_stem="$2"
output_stem="$3"
surge_results_suffix="$4"
num_jobs="$5"


# Input data
qtl_expression_file=$eqtl_input_dir"cross_tissue_eqtl_expression_input.txt"
qtl_genotype_file=$eqtl_input_dir"cross_tissue_eqtl_genotype_input.txt"
qtl_covariate_file=$eqtl_input_dir"cross_tissue_eqtl_covariate_input.txt"
qtl_test_names_file=$eqtl_input_dir"all_tests.txt"
qtl_sample_overlap_file=$eqtl_input_dir"individual_id.txt"

#######################################
threshold="1e-5"
surge_latent_factor_file=$output_stem"perm_False_surge_latent_factors_"$threshold".txt"
perm_surge_latent_factor_factor_file=$output_stem"perm_interaction_only_surge_latent_factors_"$threshold".txt"
if false; then
python extract_surge_latent_factor_by_pve.py $surge_results_stem $surge_latent_factor_file $perm_surge_latent_factor_factor_file $surge_results_suffix $threshold
fi

# Important: assuming first couple of columns of qtl_covariate_file are expression PCs
expression_pc_factor_file=$output_stem"expression_pc_factors_"$threshold".txt"
qtl_covariate_remove_top_expression_pc_file=$output_stem"cross_tissue_eqtl_covariate_without_top_expression_pcs_"$threshold"_input.txt"
if false; then
python2 extract_expression_pc_factor_file_with_matched_number_of_components.py $surge_latent_factor_file $qtl_covariate_file $expression_pc_factor_file $qtl_covariate_remove_top_expression_pc_file
fi


# Interaction analysis on real data
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"perm_False_"$threshold"_interaction_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $surge_latent_factor_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
done
fi


# Interaction analysis on permuted data
if false; then
sample_permutation_file=$surge_results_stem"interaction_only"$surge_results_suffix"var_param_1e-3_sample_permutation.txt"
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"perm_interaction_only_"$threshold"_interaction_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_permuted_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $perm_surge_latent_factor_factor_file $qtl_sample_overlap_file $qtl_output_root $sample_permutation_file $job_number $num_jobs
done
fi


# Interaction eqtl analysis on expression pc data
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"perm_False_expression_pc_"$threshold"_interaction_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_remove_top_expression_pc_file $expression_pc_factor_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
done
fi

if false; then
qtl_output_root=$output_stem"perm_False_"$threshold"_interaction_eqtl_results_"
python2 merge_parallelized_latent_factor_interaction_eqtl_calls.py $qtl_output_root $num_jobs

qtl_output_root=$output_stem"perm_interaction_only_"$threshold"_interaction_eqtl_results_"
python2 merge_parallelized_latent_factor_interaction_eqtl_calls.py $qtl_output_root $num_jobs

qtl_output_root=$output_stem"perm_False_expression_pc_"$threshold"_interaction_eqtl_results_"
python2 merge_parallelized_latent_factor_interaction_eqtl_calls.py $qtl_output_root $num_jobs
fi

#######################################
threshold="5e-4"
surge_latent_factor_file=$output_stem"perm_False_surge_latent_factors_"$threshold".txt"
perm_surge_latent_factor_factor_file=$output_stem"perm_interaction_only_surge_latent_factors_"$threshold".txt"
if false; then
python extract_surge_latent_factor_by_pve.py $surge_results_stem $surge_latent_factor_file $perm_surge_latent_factor_factor_file $surge_results_suffix $threshold
fi

# Important: assuming first couple of columns of qtl_covariate_file are expression PCs
expression_pc_factor_file=$output_stem"expression_pc_factors_"$threshold".txt"
qtl_covariate_remove_top_expression_pc_file=$output_stem"cross_tissue_eqtl_covariate_without_top_expression_pcs_"$threshold"_input.txt"
if false; then
python2 extract_expression_pc_factor_file_with_matched_number_of_components.py $surge_latent_factor_file $qtl_covariate_file $expression_pc_factor_file $qtl_covariate_remove_top_expression_pc_file
fi


# Interaction analysis on real data
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"perm_False_"$threshold"_interaction_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $surge_latent_factor_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
done
fi



# Interaction analysis on permuted data
if false; then
sample_permutation_file=$surge_results_stem"interaction_only"$surge_results_suffix"var_param_1e-3_sample_permutation.txt"
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"perm_interaction_only_"$threshold"_interaction_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_permuted_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $perm_surge_latent_factor_factor_file $qtl_sample_overlap_file $qtl_output_root $sample_permutation_file $job_number $num_jobs
done
fi


# Interaction eqtl analysis on expression pc data
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_stem"perm_False_expression_pc_"$threshold"_interaction_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_remove_top_expression_pc_file $expression_pc_factor_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
done
fi

if false; then
qtl_output_root=$output_stem"perm_False_"$threshold"_interaction_eqtl_results_"
python2 merge_parallelized_latent_factor_interaction_eqtl_calls.py $qtl_output_root $num_jobs

qtl_output_root=$output_stem"perm_interaction_only_"$threshold"_interaction_eqtl_results_"
python2 merge_parallelized_latent_factor_interaction_eqtl_calls.py $qtl_output_root $num_jobs

qtl_output_root=$output_stem"perm_False_expression_pc_"$threshold"_interaction_eqtl_results_"
python2 merge_parallelized_latent_factor_interaction_eqtl_calls.py $qtl_output_root $num_jobs
fi







module load r/3.6.3
Rscript visualize_surge_interaction_eqtls.R $output_stem"perm_"




