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
# Output root
job_number="0"
output_root=$standard_eqtl_results_dir"sc_standard_eqtl_analysis_scran_1000_hvg_"
sh run_eqtl_analysis_with_random_effects_in_parallel.sh $expression_file $genotype_file $test_names_file $covariate_file $sample_overlap_file $output_root $job_number $total_jobs


if false; then
output_root=$per_time_step_eqtl_dir"sc_per_time_step_eqtl_analysis_"$day_num"_day_"$num_pcs"_pcs_"
for job_number in $(seq 0 `expr $total_jobs - "1"`); do
	sbatch run_eqtl_analysis_with_random_effects_in_parallel.sh $expression_file $genotype_file $test_names_file $covariate_file $sample_overlap_file $num_pcs $output_root $job_number $total_jobs
done
fi
