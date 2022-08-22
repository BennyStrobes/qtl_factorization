#!/bin/bash -l
#SBATCH
#SBATCH --time=40:00:00
#SBATCH --mem=20GB
#SBATCH --nodes=1

qtl_expression_file="$1"
qtl_genotype_file="$2"
qtl_test_names_file="$3"
qtl_covariate_file="$4"
qtl_sample_overlap_file="$5"
qtl_output_root="$6"
job_number="$7"
num_jobs="$8"
version="$9"

total_lines=`wc -l $qtl_test_names_file`



module load r/3.6.3

Rscript run_eqtl_analysis_in_parallel.R $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $qtl_sample_overlap_file $qtl_output_root $version $job_number $num_jobs $total_lines

