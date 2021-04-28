#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=lrgmem
#SBATCH --mem=20GB
#SBATCH --nodes=1

tissues_file="$1"
gtex_expression_dir="$2"
gtex_tpm_dir="$3"
gtex_covariate_dir="$4"
gtex_genotype_dir="$5"
gtex_egene_dir="$6"
gtex_individual_information_file="$7"
gtex_sample_information_file="$8"
gtex_eqtl_dir="$9"
gtex_xcell_enrichment_file="${10}"
output_dir="${11}"
visualization_expression_dir="${12}"
eqtl_visualization_dir="${13}"



module load python/2.7-anaconda
##############################################
# Preprocess gene expression data
##############################################
if false; then
python preprocess_gtex_expression_data.py $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $gtex_individual_information_file $gtex_sample_information_file $gtex_eqtl_dir $gtex_xcell_enrichment_file $output_dir
fi

if false; then
module load R/3.5.1
Rscript visualize_processed_gtex_expression.R $tissues_file $output_dir $visualization_expression_dir
fi

#############################################
# Prepare data for standard eqtl analysis
#############################################
if false; then
python prepare_gtex_data_for_cross_tissue_eqtl_analysis.py $tissues_file $gtex_genotype_dir $output_dir
fi




#############################################
# Run latent-factor interaction QTL Analysis
#############################################
# Input data
# Input data
qtl_expression_file=$output_dir"cross_tissue_eqtl_residual_expression_input.txt"
qtl_genotype_file=$output_dir"cross_tissue_eqtl_genotype_input.txt"
qtl_covariate_file=$output_dir"cross_tissue_eqtl_residual_covariate_input.txt"
qtl_interaction_factor_file=$output_dir"cross_tissue_eqtl_residual_interaction_factor_input.txt"
qtl_test_names_file=$output_dir"all_tests.txt"
qtl_sample_overlap_file=$output_dir"individual_id.txt"

#In parallel
num_jobs="30"
job_number="0"
qtl_output_root=$output_dir"cross_tissue_latent_factor_interaction_eqtl_results_"$job_number"_"$num_jobs"_"
if false; then
sbatch run_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $qtl_interaction_factor_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
fi
if false; then
for job_number in $(seq 1 $(($num_jobs-1))); do 
	qtl_output_root=$output_dir"cross_tissue_latent_factor_interaction_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $qtl_interaction_factor_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
done
fi


# merge parallel runs
if false; then
python merge_parallelized_latent_factor_interaction_eqtl_calls.py $output_dir"cross_tissue_latent_factor_interaction_eqtl_results_" $num_jobs
fi

#############################################
# Run known tissue identity interaction QTL Analysis
#############################################
# Input data
qtl_expression_file=$output_dir"cross_tissue_eqtl_residual_expression_input.txt"
qtl_genotype_file=$output_dir"cross_tissue_eqtl_genotype_input.txt"
qtl_covariate_file=$output_dir"cross_tissue_eqtl_residual_covariate_known_tissue_input.txt"
qtl_interaction_factor_file=$output_dir"cross_tissue_eqtl_residual_known_tissue_interaction_input.txt"
qtl_test_names_file=$output_dir"all_tests.txt"
qtl_sample_overlap_file=$output_dir"individual_id.txt"

#In parallel
num_jobs="30"
job_number="0"
qtl_output_root=$output_dir"cross_tissue_known_tissue_interaction_eqtl_results_"$job_number"_"$num_jobs"_"
if false; then
sbatch run_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $qtl_interaction_factor_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
fi

if false; then
for job_number in $(seq 1 $(($num_jobs-1))); do 
	qtl_output_root=$output_dir"cross_tissue_known_tissue_interaction_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $qtl_interaction_factor_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
done
fi

# merge parallel runs
if false; then
python merge_parallelized_latent_factor_interaction_eqtl_calls.py $output_dir"cross_tissue_known_tissue_interaction_eqtl_results_" $num_jobs
fi


if false; then
module load R/3.5.1
Rscript visualize_standard_eqtls.R $output_dir $tissues_file $eqtl_visualization_dir 
fi

#############################################
# Prepare data for eQTL factorization 
#############################################
num_genes="2000"
python preprocess_gtex_data_based_on_latent_factor_interaction_eqtls_for_eqtl_factorization.py $output_dir $num_genes





































#########################
# OLD 
#########################


#############################################
# Run standard eqtl analysis
#############################################
# Input data
qtl_expression_file=$output_dir"cross_tissue_eqtl_expression_input.txt"
qtl_genotype_file=$output_dir"cross_tissue_eqtl_genotype_input.txt"
qtl_covariate_file=$output_dir"cross_tissue_eqtl_covariate_input.txt"
qtl_test_names_file=$output_dir"all_tests.txt"
qtl_sample_overlap_file=$output_dir"individual_id.txt"
# Output root
qtl_output_root=$output_dir"cross_tissue_eqtl_results_"
if false; then
# No random effects term
python run_eqtl_analysis.py $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $qtl_output_root
fi


# Random effects in parallel
num_jobs="75"

job_number="0"
if false; then
qtl_output_root=$output_dir"cross_tissue_re_eqtl_results_"$job_number"_"$num_jobs"_"
sbatch run_eqtl_analysis_with_random_effects_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
fi
if false; then
for job_number in $(seq 1 $(($num_jobs-1))); do 
	qtl_output_root=$output_dir"cross_tissue_re_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_eqtl_analysis_with_random_effects_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
done
fi


# merge parallel runs
if false; then
python merge_parallelized_eqtl_calls.py $output_dir"cross_tissue_re_eqtl_results_" $num_jobs
fi



num_genes="2000"
random_seed="1"
nominal_pvalue_thresh=".01"
if false; then 
python preprocess_gtex_data_for_eqtl_factorization.py $output_dir $num_genes $random_seed $nominal_pvalue_thresh
fi

num_genes="2000"
random_seed="1"
nominal_pvalue_thresh=".000000001"
if false; then
python preprocess_gtex_data_for_eqtl_factorization.py $output_dir $num_genes $random_seed $nominal_pvalue_thresh
fi


num_genes="2000"
random_seed="1"
if false; then
python preprocess_gtex_data_based_on_tissue_egenes_for_eqtl_factorization.py $output_dir $num_genes $random_seed $tissues_file $gtex_egene_dir
fi



