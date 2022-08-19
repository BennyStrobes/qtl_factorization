#!/bin/bash -l

#SBATCH
#SBATCH --time=50:00:00
#SBATCH --partition=shared
#SBATCH --mem=20GB
#SBATCH --nodes=1


#########################
# Output directories
#########################
# Root directory for simulation analysis
root_dir="/scratch16/abattle4/bstrober/qtl_factorization/gtex_cis_eqtl/"
# Directory containing input data
input_data_dir=$root_dir"input_data/"
# Directory containing simulated data
processed_data_dir=$root_dir"processed_data/"
# Directory containing eqtl results on simulated data
eqtl_results_dir=$root_dir"eqtl_results/"
# Directory containing interaction eqtl results
interaction_eqtl_dir=$root_dir"interaction_eqtls/"
# Directory containing visualizations
visualization_expression_dir=$root_dir"visualize_expression/"
# Directory containing visualizations of standard eqtl approachees
visualize_standard_eqtl_dir=$root_dir"visualize_standard_eqtl/"
# Directory containing visualization of eQTL Factorization results
visualize_eqtl_factorization_results_dir=$root_dir"visualize_eqtl_factorization/"

#########################
# Input Data
#########################
gtex_expression_dir="/scratch16/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_expression_matrices/"
gtex_tpm_dir="/scratch16/abattle4/lab_data/GTEx_v8/processed/rna_seq_by_tissue/gene_tpm/"
gtex_covariate_dir="/scratch16/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_covariates/"
gtex_genotype_dir="/scratch16/abattle4/lab_data/GTEx_v8_trans_eqtl_data_processed_by_brian/processed_genotypes/"
gtex_egene_dir="/scratch16/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL/"
gtex_eqtl_dir="/scratch16/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_all_associations/"
gtex_tissue_colors_file="/scratch16/abattle4/bstrober/qtl_factorization/gtex_cis_eqtl/input_data/gtex_tissue_colors.txt"
gtex_individual_information_file="/scratch16/abattle4/lab_data/GTEx_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"
gtex_sample_information_file="/scratch16/abattle4/bstrober/qtl_factorization/gtex_cis_eqtl/input_data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
gtex_xcell_enrichment_file="/scratch16/abattle4/bstrober/qtl_factorization/gtex_cis_eqtl/input_data/GTEx_Analysis_v8_xCell_scores_7_celltypes.txt"





#####################
## 10 tissues case
#####################
tissues_file=$input_data_dir"tissues_subset_10.txt"
output_dir=$processed_data_dir"tissues_subset_10_"
output_visualization_dir=$visualization_expression_dir"tissues_subset_10_"
output_eqtl_visualization_dir=$visualize_standard_eqtl_dir"tissues_subset_10_"
qtl_model_version="lm"
num_qtl_jobs="20"
sh preprocess_gtex_data_for_eqtl_factorization.sh $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $gtex_individual_information_file $gtex_sample_information_file $gtex_eqtl_dir $gtex_xcell_enrichment_file $output_dir $output_visualization_dir $output_eqtl_visualization_dir $qtl_model_version $num_qtl_jobs













#############################################
## Run SURGE!
#############################################
input_data_stem="tissues_subset_10_"
sample_overlap_file=$processed_data_dir$input_data_stem"individual_id.txt"
expression_training_file=$processed_data_dir$input_data_stem"eqtl_factorization_standard_eqtl_2000_lm_input_expression.npy"
genotype_training_file=$processed_data_dir$input_data_stem"eqtl_factorization_standard_eqtl_2000_lm_input_genotype.npy"
covariate_file=$processed_data_dir$input_data_stem"cross_tissue_eqtl_2000_lm_covariate_w_intercept_input.txt"
num_latent_factors="20"
lambda_v="1"
variance_param="1e-3"
ard_variance_param="1e-3"
seed="1"
model_name="surge"
ratio_variance_standardization="True"

if false; then
num_latent_factors="20"
permutation_type="False"
warmup_iterations="5"
re="True"
ard_variance_param="1e-3"
variance_param="1e-3"
seed="1"
output_stem=$eqtl_results_dir$input_data_stem$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_re_"$re"_var_param_"$ard_variance_param"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $re

num_latent_factors="20"
permutation_type="interaction_only"
warmup_iterations="5"
re="True"
ard_variance_param="1e-3"
variance_param="1e-3"
seed="1"
output_stem=$eqtl_results_dir$input_data_stem$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_re_"$re"_var_param_"$ard_variance_param"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $re

num_latent_factors="20"
permutation_type="False"
warmup_iterations="5"
re="False"
ard_variance_param="1e-3"
variance_param="1e-3"
seed="1"
output_stem=$eqtl_results_dir$input_data_stem$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_re_"$re"_var_param_"$ard_variance_param"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $re

num_latent_factors="20"
permutation_type="interaction_only"
warmup_iterations="5"
re="False"
ard_variance_param="1e-3"
variance_param="1e-3"
seed="1"
output_stem=$eqtl_results_dir$input_data_stem$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_re_"$re"_var_param_"$ard_variance_param"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $re



# Also run forvariance parameter prior of 1e-16
num_latent_factors="20"
permutation_type="False"
warmup_iterations="5"
re="True"
ard_variance_param="1e-16"
variance_param="1e-16"
seed="1"
output_stem=$eqtl_results_dir$input_data_stem$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_re_"$re"_var_param_"$ard_variance_param"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $re

num_latent_factors="20"
permutation_type="interaction_only"
warmup_iterations="5"
re="True"
ard_variance_param="1e-16"
variance_param="1e-16"
seed="1"
output_stem=$eqtl_results_dir$input_data_stem$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_re_"$re"_var_param_"$ard_variance_param"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $re

num_latent_factors="20"
permutation_type="False"
warmup_iterations="5"
re="False"
ard_variance_param="1e-16"
variance_param="1e-16"
seed="1"
output_stem=$eqtl_results_dir$input_data_stem$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_re_"$re"_var_param_"$ard_variance_param"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $re

num_latent_factors="20"
permutation_type="interaction_only"
warmup_iterations="5"
re="False"
ard_variance_param="1e-16"
variance_param="1e-16"
seed="1"
output_stem=$eqtl_results_dir$input_data_stem$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_re_"$re"_var_param_"$ard_variance_param"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $re
fi










#############################################
## Run interaction eqtl analysis with SURGE-latent factors
#############################################

# SURGE Model parameters
input_data_stem="tissues_subset_10_"
lambda_v="1"
model_name="surge"
ratio_variance_standardization="True"
num_latent_factors="20"
permutation_type="False"
warmup_iterations="5"
re="False"
ard_variance_param="1e-3"
variance_param="1e-3"
seed="1"

# eqtl input dir
eqtl_input_dir=$processed_data_dir$input_data_stem

# Surge output stem
surge_results_stem=$eqtl_results_dir$input_data_stem$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"
surge_results_suffix="_re_False_"

# How much parallelization
num_parallel_jobs="10"

# Output root
output_stem=$interaction_eqtl_dir$input_data_stem$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_"
if false; then
sh run_interaction_eqtl_analysis_with_surge_factors.sh $eqtl_input_dir $surge_results_stem $output_stem $surge_results_suffix $num_parallel_jobs
fi



#############################################
## Visualize results
#############################################
if false; then
module load R/3.5.1
Rscript visualize_eqtl_factorization.R $processed_data_dir $eqtl_results_dir $visualize_eqtl_factorization_results_dir $gtex_tissue_colors_file
fi



#####################
## Preprocess data in a single tissue
#####################
tissues_file=$input_data_dir"tissues_subset_colon_transverse.txt"
output_dir=$processed_data_dir"tissues_subset_colon_transverse_"
output_visualization_dir=$visualization_expression_dir"tissues_subset_colon_transverse_"
output_eqtl_visualization_dir=$visualize_standard_eqtl_dir"tissues_subset_colon_transverse_"
qtl_model_version="lm"
num_qtl_jobs="10"
if false; then
sh preprocess_gtex_data_for_eqtl_factorization.sh $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $gtex_individual_information_file $gtex_sample_information_file $gtex_eqtl_dir $gtex_xcell_enrichment_file $output_dir $output_visualization_dir $output_eqtl_visualization_dir $qtl_model_version $num_qtl_jobs
fi

#############################################
## Run interaction eqtl analysis with cell-type proportions
#############################################
input_data_stem="tissues_subset_colon_transverse_"
eqtl_input_dir=$processed_data_dir$input_data_stem
output_stem=$interaction_eqtl_dir"xcell_interaction_"$input_data_stem
if false; then
sh run_interaction_eqtl_analysis_with_cell_type_proportions.sh $eqtl_input_dir $output_stem
fi







#############################################
## Run SURGE in a single tissue
#############################################
input_data_stem="tissues_subset_colon_transverse_"
sample_overlap_file=$processed_data_dir$input_data_stem"individual_id.txt"
expression_training_file=$processed_data_dir$input_data_stem"eqtl_factorization_standard_eqtl_2000_input_expression.npy"
genotype_training_file=$processed_data_dir$input_data_stem"eqtl_factorization_standard_eqtl_2000_input_genotype.npy"
covariate_file=$processed_data_dir$input_data_stem"cross_tissue_eqtl_2000_covariate_w_intercept_input.txt"
num_latent_factors="10"
lambda_v="1"
variance_param="1e-3"
ard_variance_param="1e-3"
seed="1"
model_name="surge"
ratio_variance_standardization="True"
re="False"

# No permutation, ard prior
if false; then
permutation_type="False"
warmup_iterations="5"
output_stem=$eqtl_results_dir$input_data_stem$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_2000_tests_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $re

# Permutation, ard prior
permutation_type="interaction_only"
warmup_iterations="5"
output_stem=$eqtl_results_dir$input_data_stem$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_2000_tests_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $re
fi





#############################################
## Run interaction eqtl analysis with SURGE-latent factors
#############################################
# SURGE Model parameters
input_data_stem="tissues_subset_colon_transverse_"
lambda_v="1"
model_name="surge"
ratio_variance_standardization="True"
num_latent_factors="10"
permutation_type="False"
warmup_iterations="5"
re="False"
ard_variance_param="1e-3"
variance_param="1e-3"
seed="1"

# eqtl input dir
eqtl_input_dir=$processed_data_dir$input_data_stem

# Surge output stem
surge_results_stem=$eqtl_results_dir$input_data_stem$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"
surge_results_suffix="_2000_tests_"

# How much parallelization
num_parallel_jobs="20"

# Output root
output_stem=$interaction_eqtl_dir$input_data_stem$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_"
if false; then
sh run_interaction_eqtl_analysis_with_surge_factors.sh $eqtl_input_dir $surge_results_stem $output_stem $surge_results_suffix $num_parallel_jobs
fi






########################################
# Visualize results of single tissue analysis
########################################

if false; then
module load R/3.5.1
Rscript visualize_single_tissue_eqtl_factorization.R $processed_data_dir $eqtl_results_dir $visualize_eqtl_factorization_results_dir $gtex_tissue_colors_file $interaction_eqtl_dir
fi

