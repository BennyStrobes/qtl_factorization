#!/bin/bash -l

#SBATCH
#SBATCH --time=2:00:00
#SBATCH --partition=express

#SBATCH --nodes=1

######################
# Input data
######################
# Input single cell expression data (emailed by Meena Subramaniam on Nov. 18, 2019)
input_h5py_file="/scratch16/abattle4/lab_data/ye_lab_lupus_data_set/expression/Lupus_study_adjusted.h5ad"
##
# Gene annotation file (hg19)
gene_annotation_file="/scratch16/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/input_data/gencode.v19.annotation.gff3"

# Directory containing genotype data
genotype_data_dir="/scratch16/abattle4/lab_data/ye_lab_lupus_data_set/imputed_genotypes/"

# File containing mapping from old genotype ids to new genotype ids
# Basically rna-seq had old genotype ids and genotype data had new genotype ids
# We will convert rna-seq to new genotype ids
genotype_id_mapping_file="/scratch16/abattle4/lab_data/ye_lab_lupus_data_set/imputed_genotypes/sample_mappings_ye_SLE.txt"

# File containing genotyped individuals
genotyped_individuals_file="/scratch16/abattle4/lab_data/ye_lab_lupus_data_set/imputed_genotypes/genotyped_individuals.txt"

# File containing isg scores for each individual
isg_score_file="/scratch16/abattle4/lab_data/ye_lab_lupus_data_set/covariate_info/isg_scores.csv"

# File containing isg scores for each cell
cell_isg_score_file="/scratch16/abattle4/lab_data/ye_lab_lupus_data_set/covariate_info/cell_isg_scores.csv"

# Gene set enrichment data
# https://github.com/princyparsana/gsea
gsea_data_dir="/scratch16/abattle4/bstrober/tools/gsea/data/"

# Directory containing coloc input data
#coloc_input_dir="/work-zfs/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/coloc/input_data/"


# Directory containing ldsc source code
ldsc_source_code_dir="/scratch16/abattle4/bstrober/tools/ldsc/"
custom_ldsc_source_code_dir="/scratch16/abattle4/bstrober/tools/custom_sldsc_multiple_independent_annotations/"

# Directory containing Alkes group processed sumstats: Found on o2 at: /n/groups/price/ldsc/sumstats_formatted_2021/
ldsc_sumstats_dir="/scratch16/abattle4/bstrober/sumstats/alkesgrp_sumstats/"

# Directory containing UKBB sumstats: found on o2 at: /n/groups/price/UKBiobank/sumstats/sumstats_for_ukb/
# Also contains bentham et al 2015 sle sumstats found here: http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003156/
ukbb_sumstats_dir="/scratch16/abattle4/bstrober/sumstats/ukbb_sumstats/"

# Directory containing sldsc input data
sldsc_input_data_dir="/scratch16/abattle4/bstrober/ldsc_input_data/"

######################
# Output directories
######################
# Root of all output directoreis
output_root="/scratch16/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/"

# Directory containing processed genotype data
processed_genotype_dir=$output_root"processed_genotype/"

# Directory containing processed single cell expression
processed_expression_dir=$output_root"processed_expression/"

# Directory containing processed pseudobulk single cell expression
processed_pseudobulk_expression_dir=$output_root"processed_pseudobulk_expression/"

# Directory containing visualizations of processed single cell expression
visualize_processed_expression_dir=$output_root"visualize_processed_expression/"

# Directory containing visualizations of processed single cell expression
visualize_processed_pseudobulk_expression_dir=$output_root"visualize_processed_pseudobulk_expression/"

# Directory containing visualizations of processed single cell expression
visualize_processed_genotype_dir=$output_root"visualize_processed_genotype/"

# Directory containing standard eqtl dir
standard_eqtl_dir=$output_root"standard_eqtl/"

# Directory containing latent factor interaction eqtls
latent_factor_interaction_eqtl_dir=$output_root"latent_factor_interaction_eqtl/"

# Directory containing visualizations of standard eqtls
visualize_latent_factor_interaction_eqtl_dir=$output_root"visualize_latent_factor_interaction_eqtl/"

# Directory containing pre-processed eqtl factorization input files
eqtl_factorization_input_dir=$output_root"eqtl_factorization_input/"

# Directory containing eqtl factorization results
eqtl_factorization_results_dir=$output_root"eqtl_factorization_results/"

# Directory containing visualizations of eqtl factorization results
eqtl_factorization_visualization_dir=$output_root"visualize_eqtl_factorization/"


# Directory containing surge interaction eqtl results
surge_interaction_eqtl_dir=$output_root"surge_interaction_eqtl_v2/"

# Directory containing gene set enrichment results
gene_set_enrichment_dir=$output_root"gene_set_enrichment/"

# Coloc dir
coloc_dir=$output_root"coloc/"

# visualize Coloc dir
visualize_coloc_dir=$output_root"visualize_coloc/"

# sldsc processed data dir
sldsc_processed_data_dir=$output_root"sldsc_processed_data/"

# Directory containing gridspace sldsc processed data
component_gridspace_sldsc_processed_data_dir=$output_root"component_gridspace_sldsc_processed_data/"

# Directory containing gridspace sldsc results
component_gridspace_sldsc_results_dir=$output_root"component_gridspace_sldsc_results/"

# Directory containing static_eqtl sldsc processed data
static_eqtl_sldsc_processed_data_dir=$output_root"static_eqtl_sldsc_processed_data/"

# Directory containing static eqtl sldsc results
static_eqtl_sldsc_results_dir=$output_root"static_eqtl_sldsc_results/"





# Directory containing results of eQTL correlation analysis
eqtl_correlation_dir=$output_root"eqtl_correlation/"

# Directory containing per-cell SLDSC processed data
per_cell_sldsc_processed_data_dir=$output_root"per_cell_sldsc_processed_data/"

# Directory containing per-cell SLDSC results
per_cell_sldsc_results_dir=$output_root"per_cell_sldsc_results/"

# Directory containing per-cell SLDSC processed data
per_cell_3_component_sldsc_processed_data_dir=$output_root"per_cell_3_component_sldsc_processed_data/"

# Directory containing per-cell SLDSC results
per_cell_3_component_sldsc_results_dir=$output_root"per_cell_3_component_sldsc_results/"

# Directory containing gridspace sldsc processed data
component_gridspace_sldsc_processed_data_dir=$output_root"component_gridspace_sldsc_processed_data/"

# Directory containing gridspace sldsc results
component_gridspace_sldsc_results_dir=$output_root"component_gridspace_sldsc_results/"



######################
# Preprocess single cell expression
######################
regress_out_batch="True"
sh preprocess_single_cell_expression.sh $input_h5py_file $processed_expression_dir $visualize_processed_expression_dir $gene_annotation_file $genotyped_individuals_file $processed_pseudobulk_expression_dir $visualize_processed_pseudobulk_expression_dir $regress_out_batch $isg_score_file $cell_isg_score_file


######################
# Preprocess Genotype data
######################
# file containing list of individuals with have sc rna-seq for (generated by preprocess_single_cell_expression.sh)
sample_covariates_file=$processed_pseudobulk_expression_dir"no_outlier_pseudobulk_scran_normalization_hvg_6000_regress_batch_True_individual_clustering_leiden_resolution_10.0_no_cap_15_sample_covariates.txt"
if false; then
sh preprocess_genotype.sh $genotype_data_dir $processed_genotype_dir $sample_covariates_file $visualize_processed_genotype_dir
fi


######################
# Call latent-factor (PCA) interaction-eqtls
######################
if false; then
sh latent_factor_interaction_eqtl_driver_key.sh $processed_expression_dir $processed_pseudobulk_expression_dir $processed_genotype_dir $gene_annotation_file $latent_factor_interaction_eqtl_dir $visualize_latent_factor_interaction_eqtl_dir
fi

#############################################
# Prepare data for eQTL factorization 
#############################################
if false; then
sh preprocess_data_for_eqtl_factorization.sh $latent_factor_interaction_eqtl_dir $eqtl_factorization_input_dir
fi


#############################################
## Run eqtl factorization!
#############################################
input_data_stem="eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore"
sample_overlap_file=$eqtl_factorization_input_dir$input_data_stem"_eqtl_input_sample_overlap.txt"
expression_training_file=$eqtl_factorization_input_dir$input_data_stem"_eqtl_input_expression.npy"
genotype_training_file=$eqtl_factorization_input_dir$input_data_stem"_eqtl_input_unnormalized_genotype.npy"
covariate_file=$eqtl_factorization_input_dir$input_data_stem"_eqtl_input_covariates.txt"
test_names_file=$eqtl_factorization_input_dir$input_data_stem"_eqtl_input_test_names.txt"
num_latent_factors="10"
lambda_v="1"
variance_param="1e-3"
ard_variance_param="1e-3"
seed="1"

if false; then

permutation_type="False"
model_name="surge"
ratio_variance_standardization="False"
round_genotype="False"
warmup_iterations="5"
num_latent_factors="10"
data_filter="filter_RP"
delta_elbo_thresh="1e-3"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_"$num_latent_factors"_seed_"$seed"_warm_"$warmup_iterations"_rv_std_"$ratio_variance_standardization"_perm_"$permutation_type"_delta_elbo_"$delta_elbo_thresh"_alt_init_"

sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype $data_filter $test_names_file $delta_elbo_thresh


permutation_type="False"
model_name="surge"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="5"
num_latent_factors="10"
data_filter="filter_RP"
delta_elbo_thresh="1e-1"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_"$num_latent_factors"_seed_"$seed"_warm_"$warmup_iterations"_rv_std_"$ratio_variance_standardization"_perm_"$permutation_type"_delta_elbo_"$delta_elbo_thresh"_alt_init_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype $data_filter $test_names_file $delta_elbo_thresh
fi
if false; then
permutation_type="False"
model_name="surge"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="5"
num_latent_factors="10"
data_filter="filter_RP"
delta_elbo_thresh="5e-1"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_"$num_latent_factors"_seed_"$seed"_warm_"$warmup_iterations"_rv_std_"$ratio_variance_standardization"_perm_"$permutation_type"_delta_elbo_"$delta_elbo_thresh"_alt_init_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype $data_filter $test_names_file $delta_elbo_thresh
fi
if false; then

permutation_type="interaction_only"
model_name="surge"
ratio_variance_standardization="False"
round_genotype="False"
warmup_iterations="5"
num_latent_factors="10"
data_filter="filter_RP"
delta_elbo_thresh="1e-3"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_"$num_latent_factors"_seed_"$seed"_warm_"$warmup_iterations"_rv_std_"$ratio_variance_standardization"_perm_"$permutation_type"_delta_elbo_"$delta_elbo_thresh"_alt_init_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype $data_filter $test_names_file $delta_elbo_thresh

permutation_type="interaction_only"
model_name="surge"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="5"
num_latent_factors="10"
data_filter="filter_RP"
delta_elbo_thresh="1e-1"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_"$num_latent_factors"_seed_"$seed"_warm_"$warmup_iterations"_rv_std_"$ratio_variance_standardization"_perm_"$permutation_type"_delta_elbo_"$delta_elbo_thresh"_alt_init_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype $data_filter $test_names_file $delta_elbo_thresh
fi
if false; then

permutation_type="interaction_only"
model_name="surge"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="5"
num_latent_factors="10"
data_filter="filter_RP"
delta_elbo_thresh="5e-1"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_"$num_latent_factors"_seed_"$seed"_warm_"$warmup_iterations"_rv_std_"$ratio_variance_standardization"_perm_"$permutation_type"_delta_elbo_"$delta_elbo_thresh"_alt_init_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype $data_filter $test_names_file $delta_elbo_thresh
fi


if false; then
permutation_type="False"
model_name="surge"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="5"
num_latent_factors="10"
data_filter="filter_hwe"
delta_elbo_thresh="1e-2"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_"$num_latent_factors"_seed_"$seed"_warm_"$warmup_iterations"_rv_std_"$ratio_variance_standardization"_perm_"$permutation_type"_delta_elbo_"$delta_elbo_thresh"_"$data_filter"_alt_init_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype $data_filter $test_names_file $delta_elbo_thresh

permutation_type="interaction_only"
model_name="surge"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="5"
num_latent_factors="10"
data_filter="filter_hwe"
delta_elbo_thresh="1e-2"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_"$num_latent_factors"_seed_"$seed"_warm_"$warmup_iterations"_rv_std_"$ratio_variance_standardization"_perm_"$permutation_type"_delta_elbo_"$delta_elbo_thresh"_"$data_filter"_alt_init_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype $data_filter $test_names_file $delta_elbo_thresh

permutation_type="False"
model_name="surge"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="5"
num_latent_factors="10"
data_filter="filter_very_common"
delta_elbo_thresh="1e-2"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_"$num_latent_factors"_seed_"$seed"_warm_"$warmup_iterations"_rv_std_"$ratio_variance_standardization"_perm_"$permutation_type"_delta_elbo_"$delta_elbo_thresh"_"$data_filter"_alt_init_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype $data_filter $test_names_file $delta_elbo_thresh

permutation_type="interaction_only"
model_name="surge"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="5"
num_latent_factors="10"
data_filter="filter_very_common"
delta_elbo_thresh="1e-2"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_"$num_latent_factors"_seed_"$seed"_warm_"$warmup_iterations"_rv_std_"$ratio_variance_standardization"_perm_"$permutation_type"_delta_elbo_"$delta_elbo_thresh"_"$data_filter"_alt_init_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype $data_filter $test_names_file $delta_elbo_thresh



permutation_type="False"
model_name="surge"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="5"
num_latent_factors="10"
data_filter="filter_hwe"
delta_elbo_thresh="1e-3"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_"$num_latent_factors"_seed_"$seed"_warm_"$warmup_iterations"_rv_std_"$ratio_variance_standardization"_perm_"$permutation_type"_delta_elbo_"$delta_elbo_thresh"_"$data_filter"_alt_init_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype $data_filter $test_names_file $delta_elbo_thresh

permutation_type="interaction_only"
model_name="surge"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="5"
num_latent_factors="10"
data_filter="filter_hwe"
delta_elbo_thresh="1e-3"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_"$num_latent_factors"_seed_"$seed"_warm_"$warmup_iterations"_rv_std_"$ratio_variance_standardization"_perm_"$permutation_type"_delta_elbo_"$delta_elbo_thresh"_"$data_filter"_alt_init_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype $data_filter $test_names_file $delta_elbo_thresh

permutation_type="False"
model_name="surge"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="5"
num_latent_factors="10"
data_filter="filter_very_common"
delta_elbo_thresh="1e-3"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_"$num_latent_factors"_seed_"$seed"_warm_"$warmup_iterations"_rv_std_"$ratio_variance_standardization"_perm_"$permutation_type"_delta_elbo_"$delta_elbo_thresh"_"$data_filter"_alt_init_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype $data_filter $test_names_file $delta_elbo_thresh

permutation_type="interaction_only"
model_name="surge"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="5"
num_latent_factors="10"
data_filter="filter_very_common"
delta_elbo_thresh="1e-3"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_"$num_latent_factors"_seed_"$seed"_warm_"$warmup_iterations"_rv_std_"$ratio_variance_standardization"_perm_"$permutation_type"_delta_elbo_"$delta_elbo_thresh"_"$data_filter"_alt_init_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype $data_filter $test_names_file $delta_elbo_thresh
fi










############################################
# Compute interaction eqtls with surge latent factors
############################################
# Output root
output_stem=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_"
output_stem=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_v2_"  # THIS CORRESPONDS TO THE PARR RUN
output_stem=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_v3_"  # THIS CORRESPONDS TO THE HWE RUN
if false; then
sh run_interaction_eqtl_analysis_with_surge_factors.sh $latent_factor_interaction_eqtl_dir"latent_factor_interaction_hvg_6000_10.0_no_cap_15_cis_window_200000_geno_filter_False_none_zscore_eqtl_input_" $eqtl_factorization_results_dir $output_stem
fi





############################################
# Run gene set enrichment analysis
############################################
model_loadings_file=$eqtl_factorization_results_dir"eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_surge_results_k_10_seed_1_warm_5_rv_std_True_perm_False_delta_elbo_1e-2_filter_hwe_alt_init_U_S.txt"
model_pve_file=$eqtl_factorization_results_dir"eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_surge_results_k_10_seed_1_warm_5_rv_std_True_perm_False_delta_elbo_1e-2_filter_hwe_alt_init_factor_pve.txt"
expr_file=$processed_pseudobulk_expression_dir"no_outlier_pseudobulk_scran_normalization_hvg_6000_regress_batch_True_individual_clustering_leiden_resolution_10.0_no_cap_15_none_sample_norm_zscore_gene_norm_normalized_expression.txt"
gene_names_file=$processed_pseudobulk_expression_dir"no_outlier_pseudobulk_scran_normalization_hvg_6000_regress_batch_True_individual_clustering_leiden_resolution_10.0_no_cap_15_gene_names.txt"
gene_set_enrichment_output_stem=$gene_set_enrichment_dir"eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_surge_results_k_10_seed_1_warm_5_rv_std_True_perm_False_delta_elbo_1e-2_filter_hwe_alt_init_"
if false; then
sh gene_set_enrichment_analysis.sh $model_loadings_file $model_pve_file $expr_file $gene_names_file $gsea_data_dir $gene_set_enrichment_output_stem
fi



############################################
# Run component-gridspace S-LDSC
############################################
loading_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_v3_perm_False_surge_latent_factors.txt"
surge_eqtl_effect_sizes_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_v3_perm_False_interaction_eqtl_results_betas_merged.txt"
if false; then
sh run_component_gridspace_sldsc_analysis.sh $ldsc_source_code_dir $sldsc_input_data_dir $sldsc_processed_data_dir $component_gridspace_sldsc_processed_data_dir $component_gridspace_sldsc_results_dir $custom_ldsc_source_code_dir $loading_file $surge_eqtl_effect_sizes_file $ukbb_sumstats_dir $ldsc_sumstats_dir
fi


############################################
# Run static eqtl S-LDSC
############################################
static_eqtl_effect_sizes_file=$latent_factor_interaction_eqtl_dir"standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_results_merged.txt"
if false; then
sh run_static_eqtl_sldsc_analysis.sh $ldsc_source_code_dir $sldsc_input_data_dir $sldsc_processed_data_dir $static_eqtl_sldsc_processed_data_dir $static_eqtl_sldsc_results_dir $static_eqtl_effect_sizes_file
fi




############################################
# Check for overlap with coloc
############################################
if false; then
sh run_coloc_analysis.sh $surge_interaction_eqtl_dir $ukbb_sumstats_dir $coloc_dir $visualize_coloc_dir
fi






#############################################
## Visualize eqtl factorization
#############################################
module load r/3.6.3
model_stem="eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_surge_results_k_10_seed_1_warm_5_rv_std_True_perm_False_delta_elbo_1e-2_filter_hwe_alt_init_"
output_stem="eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_surge_results_k_10_seed_1_warm_5_rv_std_True_perm_False_delta_elbo_1e-2_filter_hwe_alt_init_"
Rscript visualize_eqtl_factorization.R $processed_pseudobulk_expression_dir $eqtl_factorization_results_dir $eqtl_factorization_visualization_dir $model_stem $output_stem $per_cell_sldsc_results_dir $per_cell_3_component_sldsc_results_dir $component_gridspace_sldsc_results_dir $static_eqtl_sldsc_results_dir



























if false; then
sh run_sldsc_analysis.sh $ldsc_source_code_dir $sldsc_input_data_dir $ukbb_sumstats_dir $sldsc_processed_data_dir $sldsc_results_dir $sldsc_visualization_dir $surge_interaction_eqtl_dir $ldsc_sumstats_dir
fi












############################################
# Correlation of eQTLs across components
############################################
if false; then
sh run_eqtl_correlation_analysis.sh $surge_interaction_eqtl_dir $eqtl_correlation_dir
fi





############################################
# Run per-cell S-LDSC
############################################
sample_names_file=$processed_pseudobulk_expression_dir"pseudobulk_scran_normalization_hvg_6000_regress_batch_True_individual_clustering_leiden_resolution_10.0_no_cap_15_sample_names.txt"
loading_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_surge_latent_factors_v2.txt"
surge_eqtl_effect_sizes_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_interaction_eqtl_standardized_genotype_results_v2_betas_merged.txt"
if false; then
sh run_per_cell_sldsc_analysis.sh $ldsc_source_code_dir $sldsc_input_data_dir $sldsc_processed_data_dir $surge_interaction_eqtl_dir $sample_names_file $per_cell_sldsc_processed_data_dir $per_cell_sldsc_results_dir $custom_ldsc_source_code_dir $loading_file $surge_eqtl_effect_sizes_file
fi


sample_names_file=$processed_pseudobulk_expression_dir"pseudobulk_scran_normalization_hvg_6000_regress_batch_True_individual_clustering_leiden_resolution_10.0_no_cap_15_sample_names.txt"
loading_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_surge_latent_factors_3_components_v2.txt"
surge_eqtl_effect_sizes_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_interaction_eqtl_standardized_genotype_results_3_surge_contexts_betas_merged.txt"
if false; then
sh run_per_cell_sldsc_analysis.sh $ldsc_source_code_dir $sldsc_input_data_dir $sldsc_processed_data_dir $surge_interaction_eqtl_dir $sample_names_file $per_cell_3_component_sldsc_processed_data_dir $per_cell_3_component_sldsc_results_dir $custom_ldsc_source_code_dir $loading_file $surge_eqtl_effect_sizes_file
fi

############################################
# Run component-gridspace S-LDSC
############################################
loading_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_surge_latent_factors_3_components_v2.txt"
surge_eqtl_effect_sizes_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_interaction_eqtl_standardized_genotype_results_3_surge_contexts_betas_merged.txt"
if false; then
sh run_component_gridspace_sldsc_analysis.sh $ldsc_source_code_dir $sldsc_input_data_dir $sldsc_processed_data_dir $component_gridspace_sldsc_processed_data_dir $component_gridspace_sldsc_results_dir $custom_ldsc_source_code_dir $loading_file $surge_eqtl_effect_sizes_file 
fi















