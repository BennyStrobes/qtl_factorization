#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


qtl_expression_file="$1"
qtl_genotype_file="$2"
qtl_test_names_file="$3"
qtl_covariate_file="$4"
qtl_sample_overlap_file="$5"
qtl_output_root="$6"
job_number="$7"
num_jobs="$8"


total_lines_plus_one=`wc -l $qtl_test_names_file`


module load R/3.5.1
Rscript run_eqtl_analysis_with_random_effects_in_parallel.R $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs $total_lines_plus_one