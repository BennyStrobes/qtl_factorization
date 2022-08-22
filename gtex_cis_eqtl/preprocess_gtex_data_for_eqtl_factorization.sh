#!/bin/bash -l
#SBATCH
#SBATCH --time=20:00:00
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
qtl_model_version="${14}"
num_jobs="${15}"


##############################################
# Preprocess gene expression data
##############################################
if false; then
python2 preprocess_gtex_expression_data.py $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $gtex_individual_information_file $gtex_sample_information_file $gtex_eqtl_dir $gtex_xcell_enrichment_file $output_dir


python2 preprocess_gtex_expression_data_with_outlier_detection.py $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $gtex_individual_information_file $gtex_sample_information_file $gtex_eqtl_dir $gtex_xcell_enrichment_file $output_dir
fi

if false; then
module load r/3.6.3
Rscript visualize_processed_gtex_expression.R $tissues_file $output_dir $visualization_expression_dir
fi


#############################################
# Prepare data for standard eqtl analysis
#############################################
if false; then
python2 prepare_gtex_data_for_cross_tissue_eqtl_analysis.py $tissues_file $gtex_genotype_dir $output_dir
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
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	qtl_output_root=$output_dir"cross_tissue_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs $qtl_model_version
done
fi


# merge parallel runs
if false; then
python2 merge_parallelized_eqtl_calls.py $output_dir"cross_tissue_eqtl_results_" $num_jobs $qtl_model_version
fi




#############################################
# Prepare data for eQTL factorization 
#############################################
num_genes="2000"
if false; then
python2 preprocess_gtex_data_based_on_standard_eqtls_for_eqtl_factorization.py $output_dir $num_genes $qtl_model_version
fi
















