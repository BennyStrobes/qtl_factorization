#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=lrgmem
#SBATCH --mem=80GB

#SBATCH --nodes=1


processed_expression_dir="$1"
processed_genotype_dir="$2"
gene_annotation_file="$3"
standard_eqtl_input_dir="$4"
standard_eqtl_results_dir="$5"


if false; then
python prepare_standard_eqtl_input.py $processed_expression_dir $processed_genotype_dir $gene_annotation_file $standard_eqtl_input_dir
fi


###########################
# Run single-cell eqtl analysis
###########################

total_jobs="20"

expression_file=$standard_eqtl_input_dir"scran_1000_hvg_eqtl_input_expression.txt"
genotype_file=$standard_eqtl_input_dir"scran_1000_hvg_eqtl_input_genotype.txt"
test_names_file=$standard_eqtl_input_dir"scran_1000_hvg_eqtl_input_variant_gene_pairs.txt"
sample_overlap_file=$standard_eqtl_input_dir"scran_1000_hvg_eqtl_input_sample_overlap.txt"
covariate_file=$standard_eqtl_input_dir"scran_1000_hvg_eqtl_input_covariates.txt"

if false; then
for job_number in $(seq 0 `expr $total_jobs - "1"`); do
	output_root=$standard_eqtl_results_dir"sc_standard_eqtl_analysis_scran_1000_hvg_eqtl_results_"$job_number"_"$total_jobs"_"
	sbatch run_eqtl_analysis_with_random_effects_in_parallel.sh $expression_file $genotype_file $test_names_file $covariate_file $sample_overlap_file $output_root $job_number $total_jobs
done
fi



###########################
# Merge results
###########################
if false; then
python merge_parallelized_standard_eqtl_calls.py $standard_eqtl_results_dir"sc_standard_eqtl_analysis_scran_1000_hvg_eqtl_results_" $total_jobs
fi



