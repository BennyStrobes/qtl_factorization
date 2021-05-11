#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --mem=250GB
#SBATCH --partition=lrgmem
#SBATCH --nodes=1


processed_expression_dir="$1"
processed_genotype_dir="$2"
gene_annotation_file="$3"
latent_factor_interaction_eqtl_dir="$4"






#############################################
# Prepare input data for latent factor interaction eqtl analysis
#############################################
# This first part needs 250GB of memory to run..
if false; then
python prepare_latent_factor_interaction_eqtl_data.py $processed_expression_dir $processed_genotype_dir $gene_annotation_file $latent_factor_interaction_eqtl_dir
fi


#############################################
# Run latent-factor interaction QTL Analysis
#############################################
# Input data
qtl_expression_file=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_eqtl_input_expression.txt"
qtl_genotype_file=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_eqtl_input_genotype.txt"
qtl_covariate_file=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_eqtl_input_covariates.txt"
qtl_interaction_factor_file=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_eqtl_input_latent_factors.txt"
qtl_test_names_file=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_eqtl_input_variant_gene_pairs.txt"
qtl_sample_overlap_file=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_eqtl_input_sample_overlap.txt"



#In parallel
num_jobs="50"
job_number="0"
qtl_output_root=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_eqtl_results_"$job_number"_"$num_jobs"_"
if false; then
sbatch run_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $qtl_interaction_factor_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
fi

if false; then
for job_number in $(seq 1 $(($num_jobs-1))); do 
	qtl_output_root=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $qtl_interaction_factor_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
done
fi







