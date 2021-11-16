#!/bin/bash -l

#SBATCH
#SBATCH --time=40:00:00
#SBATCH --partition=shared


qtl_expression_file="$1"
qtl_genotype_file="$2"
qtl_test_names_file="$3"
qtl_covariate_file="$4"
qtl_interaction_factor_file="$5"
qtl_sample_overlap_file="$6"
sample_permutation_file="$7"
qtl_output_root="$8"



module load R/3.5.1
Rscript run_permuted_latent_factor_interaction_eqtl_analysis_in_parallel.R $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $qtl_interaction_factor_file $qtl_sample_overlap_file $sample_permutation_file $qtl_output_root