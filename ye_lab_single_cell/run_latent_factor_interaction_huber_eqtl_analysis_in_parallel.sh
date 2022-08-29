#!/bin/bash -l

#SBATCH
#SBATCH --time=50:00:00
#SBATCH --nodes=1
#SBATCH --mem=20GB


qtl_expression_file="$1"
qtl_genotype_file="$2"
qtl_test_names_file="$3"
qtl_covariate_file="$4"
qtl_interaction_factor_file="$5"
qtl_sample_overlap_file="$6"
qtl_output_root="$7"
job_number="$8"
num_jobs="${9}"

module load r/3.6.3

total_lines=`wc -l $qtl_test_names_file`


Rscript run_latent_factor_interaction_huber_eqtl_analysis_in_parallel.R $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $qtl_interaction_factor_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs $total_lines
