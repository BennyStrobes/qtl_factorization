#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --mem=100GB
#SBATCH --nodes=1


processed_expression_dir="$1"
processed_pseudobulk_expression_dir="$2"
processed_genotype_dir="$3"
gene_annotation_file="$4"
latent_factor_interaction_eqtl_dir="$5"
visualize_latent_factor_interaction_eqtl_dir="$6"






#############################################
# Prepare input data for latent factor interaction eqtl analysis
#############################################
# This first part needs 250GB of memory to run..
if false; then
python2 prepare_latent_factor_interaction_eqtl_data.py $processed_pseudobulk_expression_dir $processed_genotype_dir $gene_annotation_file $latent_factor_interaction_eqtl_dir
fi

#############################################
# Run standard eQTL Analysis
#############################################
# Input data
qtl_expression_file=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_hvg_6000_10.0_no_cap_15_cis_window_200000_geno_filter_False_none_zscore_eqtl_input_expression.txt"
qtl_genotype_file=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_hvg_6000_10.0_no_cap_15_cis_window_200000_geno_filter_False_none_zscore_eqtl_input_genotype.txt"
qtl_covariate_file=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_hvg_6000_10.0_no_cap_15_cis_window_200000_geno_filter_False_none_zscore_eqtl_input_all_covariates.txt"
qtl_test_names_file=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_hvg_6000_10.0_no_cap_15_cis_window_200000_geno_filter_False_none_zscore_eqtl_input_variant_gene_pairs.txt"
qtl_sample_overlap_file=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_hvg_6000_10.0_no_cap_15_cis_window_200000_geno_filter_False_none_zscore_eqtl_input_sample_overlap.txt"



#In parallel
num_jobs="100"
job_number="0"
version="lmm"
qtl_output_root=$latent_factor_interaction_eqtl_dir"standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_results_"$job_number"_"$num_jobs"_"
if false; then
sbatch run_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs $version


for job_number in $(seq 1 $(($num_jobs-1))); do 
	qtl_output_root=$latent_factor_interaction_eqtl_dir"standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs $version
done
fi


# merge parallel runs
if false; then
python merge_parallelized_standard_eqtl_calls.py $latent_factor_interaction_eqtl_dir"standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_results_" $num_jobs
fi


#############################################
# Visualize latent factor interaction results
#############################################
if false; then
module load R/3.5.1
Rscript visualize_latent_factor_interaction_eqtls.R $latent_factor_interaction_eqtl_dir $visualize_latent_factor_interaction_eqtl_dir 
fi






























#############################################
# Run latent-factor interaction QTL Analysis
#############################################
# Input data
qtl_expression_file=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_10.0_none_zscore_eqtl_input_expression.txt"
qtl_genotype_file=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_10.0_none_zscore_eqtl_input_genotype.txt"
qtl_covariate_file=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_10.0_none_zscore_eqtl_input_covariates.txt"
qtl_interaction_factor_file=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_10.0_none_zscore_eqtl_input_latent_factors.txt"
qtl_test_names_file=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_10.0_none_zscore_eqtl_input_variant_gene_pairs.txt"
qtl_sample_overlap_file=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_10.0_none_zscore_eqtl_input_sample_overlap.txt"



#In parallel
num_jobs="50"
job_number="0"
qtl_output_root=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_10.0_none_zscore_eqtl_results_"$job_number"_"$num_jobs"_"
if false; then
sbatch run_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $qtl_interaction_factor_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
fi

if false; then
for job_number in $(seq 1 $(($num_jobs-1))); do 
	qtl_output_root=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_10.0_none_zscore_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $qtl_interaction_factor_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
done
fi


# merge parallel runs
if false; then
python merge_parallelized_latent_factor_interaction_eqtl_calls.py $latent_factor_interaction_eqtl_dir"latent_factor_interaction_10.0_none_zscore_eqtl_results_" $num_jobs
fi




