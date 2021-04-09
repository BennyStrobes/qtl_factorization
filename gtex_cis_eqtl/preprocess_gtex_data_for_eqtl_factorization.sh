#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
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
gtex_eqtl_dir="$8"
output_dir="$9"


module load python/2.7-anaconda
##############################################
# Preprocess gene expression data
##############################################
if false; then
python preprocess_gtex_expression_data.py $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $gtex_individual_information_file $gtex_eqtl_dir $output_dir
fi


#############################################
# Prepare data for standard eqtl analysis
#############################################
if false; then
python prepare_gtex_data_for_cross_tissue_eqtl_analysis.py $tissues_file $gtex_genotype_dir $output_dir
fi



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



#############################################
# Prepare data for eQTL factorization 
#############################################
num_genes="2000"
random_seed="1"
nominal_pvalue_thresh=".01"
python preprocess_gtex_data_for_eqtl_factorization.py $output_dir $num_genes $random_seed $nominal_pvalue_thresh





