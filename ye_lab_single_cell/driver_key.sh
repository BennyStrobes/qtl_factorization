#!/bin/bash -l

#SBATCH
#SBATCH --time=1:00:00
#SBATCH --partition=shared
#SBATCH --mem=30GB

#SBATCH --nodes=1

######################
# Input data
######################
# Input single cell expression data (emailed by Meena Subramaniam on Nov. 18, 2019)
input_h5py_file="/work-zfs/abattle4/lab_data/ye_lab_lupus_data_set/expression/Lupus_study_adjusted.h5ad"

# Gene annotation file (hg19)
gene_annotation_file="/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt"

# Directory containing genotype data
genotype_data_dir="/work-zfs/abattle4/lab_data/ye_lab_lupus_data_set/imputed_genotypes/"

# File containing mapping from old genotype ids to new genotype ids
# Basically rna-seq had old genotype ids and genotype data had new genotype ids
# We will convert rna-seq to new genotype ids
genotype_id_mapping_file="/work-zfs/abattle4/lab_data/ye_lab_lupus_data_set/imputed_genotypes/sample_mappings_ye_SLE.txt"

# File containing genotyped individuals
genotyped_individuals_file="/work-zfs/abattle4/lab_data/ye_lab_lupus_data_set/imputed_genotypes/genotyped_individuals.txt"

# File containing isg scores for each individual
isg_score_file="/work-zfs/abattle4/lab_data/ye_lab_lupus_data_set/covariate_info/isg_scores.csv"

# File containing isg scores for each cell
cell_isg_score_file="/work-zfs/abattle4/lab_data/ye_lab_lupus_data_set/covariate_info/cell_isg_scores.csv"

# GSEA data
gsea_data_dir="/work-zfs/abattle4/bstrober/tools/tools-master/gsea/data/"

# File containing list of GWAS files (downloaded here: https://alkesgroup.broadinstitute.org/sumstats_formatted/ on 9/15/21)
gwas_files="/work-zfs/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/input_data/ukbb_data/gwas_files.txt"

# directory containing gwas snp info
gwas_snp_info_dir="/work-zfs/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/input_data/ukbb_genotype/"

# Directory containing coloc input data
coloc_input_dir="/work-zfs/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/coloc/input_data/"

# File containing remap 2022 tfbs
# downloaded here https://remap.univ-amu.fr/download_page on 9/27/21
remap_tfbs_file="/work-zfs/abattle4/lab_data/tfbs_remap_2022/remap2022_nr_macs2_hg19_v1_0.bed"

######################
# Output directories
######################
# Root of all output directoreis
output_root="/work-zfs/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/"

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

# Directory containing gene set enrichment results
gene_set_enrichment_dir=$output_root"gene_set_enrichment/"

# Directory containing surge interaction eqtl results
surge_interaction_eqtl_dir=$output_root"surge_interaction_eqtl_v2/"

# Gwas overlap dir
gwas_overlap_dir=$output_root"gwas_overlap/"

# Coloc dir
coloc_dir=$output_root"coloc/"

# Susie input dir
susie_input_data_dir=$output_root"susie_input/"

# susie results dir
susie_results_dir=$output_root"susie_results/"

# susie visualization dir
susie_visualization_dir=$output_root"susie_visualization/"

# remap enrichment dir
remap_tfbs_dir=$output_root"remap_tfbs/"

# visualize Coloc dir
visualize_coloc_dir=$output_root"visualize_coloc/"



######################
# Preprocess single cell expression
######################
regress_out_batch="True"
if false; then
sh preprocess_single_cell_expression.sh $input_h5py_file $processed_expression_dir $visualize_processed_expression_dir $gene_annotation_file $genotyped_individuals_file $processed_pseudobulk_expression_dir $visualize_processed_pseudobulk_expression_dir $regress_out_batch $isg_score_file $cell_isg_score_file
fi

######################
# Preprocess Genotype data
######################
# file containing list of individuals with have sc rna-seq for (generated by preprocess_single_cell_expression.sh)
sample_covariates_file=$processed_pseudobulk_expression_dir"pseudobulk_scran_normalization_hvg_6000_regress_batch_True_individual_clustering_leiden_resolution_10.0_no_cap_15_sample_covariates.txt"
if false; then
sh preprocess_genotype.sh $genotype_data_dir $processed_genotype_dir $sample_covariates_file $visualize_processed_genotype_dir
fi


######################
# Call latent-factor (PCA) interaction-eqtls
######################
if false; then
sbatch latent_factor_interaction_eqtl_driver_key.sh $processed_expression_dir $processed_pseudobulk_expression_dir $processed_genotype_dir $gene_annotation_file $latent_factor_interaction_eqtl_dir $visualize_latent_factor_interaction_eqtl_dir
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
input_data_stem="eqtl_factorization_standard_eqtl_hvg_6000_10.0_none_zscore_capped"
input_data_stem="eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore"
sample_overlap_file=$eqtl_factorization_input_dir$input_data_stem"_eqtl_input_sample_overlap.txt"
expression_training_file=$eqtl_factorization_input_dir$input_data_stem"_eqtl_input_expression.npy"
genotype_training_file=$eqtl_factorization_input_dir$input_data_stem"_eqtl_input_unnormalized_genotype.npy"
covariate_file=$eqtl_factorization_input_dir$input_data_stem"_eqtl_input_covariates.txt"
num_latent_factors="10"
lambda_v="1"
variance_param="1e-3"
ard_variance_param="1e-16"
seed="1"

if false; then
permutation_type="False"
model_name="eqtl_factorization_vi_ard"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="3000"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_lambda_"$lambda_v"_round_geno_"$round_genotype"_alt_init_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype



permutation_type="False"
model_name="eqtl_factorization_vi_ard_full_component_update"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="3000"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_lambda_"$lambda_v"_round_geno_"$round_genotype"_alt_init_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype
fi


if false; then
permutation_type="False"
model_name="eqtl_factorization_vi_ard"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="0"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_lambda_"$lambda_v"_round_geno_"$round_genotype"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype

permutation_type="False"
model_name="eqtl_factorization_vi_ard"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="5"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_lambda_"$lambda_v"_round_geno_"$round_genotype"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype


permutation_type="False"
model_name="eqtl_factorization_vi_ard"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="3000"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_lambda_"$lambda_v"_round_geno_"$round_genotype"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype



permutation_type="False"
model_name="eqtl_factorization_vi_ard_full_component_update"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="0"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_lambda_"$lambda_v"_round_geno_"$round_genotype"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype

permutation_type="False"
model_name="eqtl_factorization_vi_ard_full_component_update"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="5"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_lambda_"$lambda_v"_round_geno_"$round_genotype"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype


permutation_type="False"
model_name="eqtl_factorization_vi_ard_full_component_update"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="3000"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_lambda_"$lambda_v"_round_geno_"$round_genotype"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype


permutation_type="False"
model_name="eqtl_factorization_vi_ard_environmental_fixed_effect"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="0"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_lambda_"$lambda_v"_round_geno_"$round_genotype"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype

permutation_type="False"
model_name="eqtl_factorization_vi_ard_environmental_fixed_effect"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="5"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_lambda_"$lambda_v"_round_geno_"$round_genotype"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype


permutation_type="False"
model_name="eqtl_factorization_vi_ard_environmental_fixed_effect"
ratio_variance_standardization="True"
round_genotype="False"
warmup_iterations="3000"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_lambda_"$lambda_v"_round_geno_"$round_genotype"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype




permutation_type="False"
model_name="eqtl_factorization_vi_ard_heteroskedastic"
ratio_variance_standardization="True"
round_genotype="True"
warmup_iterations="5"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_lambda_"$lambda_v"_round_geno_"$round_genotype"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype
fi
permutation_type="False"
model_name="eqtl_factorization_vi_ard_heteroskedastic"
ratio_variance_standardization="True"
round_genotype="True"
warmup_iterations="3000"
num_latent_factors="30"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_lambda_"$lambda_v"_round_geno_"$round_genotype"_"
if false; then
sh run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype
fi


permutation_type="False"
model_name="factorization_vi_ard"
ratio_variance_standardization="True"
round_genotype="True"
warmup_iterations="3000"
num_latent_factors="30"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_lambda_"$lambda_v"_round_geno_"$round_genotype"_"
if false; then
sh run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype
fi






############################################
# Compute interaction eqtls with surge latent factors
############################################
# Surge output files
surge_latent_factors_file=$eqtl_factorization_results_dir"eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_factorization_vi_ard_heteroskedastic_results_k_init_10_seed_1_warmup_3000_ratio_variance_std_True_permute_False_lambda_1_round_geno_True_temper_U_S.txt"
factor_pve_file=$eqtl_factorization_results_dir"eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_factorization_vi_ard_heteroskedastic_results_k_init_10_seed_1_warmup_3000_ratio_variance_std_True_permute_False_lambda_1_round_geno_True_temper_factor_pve.txt"
# Output root
output_stem=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_"
if false; then
sh run_interaction_eqtl_analysis_with_surge_factors.sh $latent_factor_interaction_eqtl_dir"latent_factor_interaction_hvg_6000_10.0_no_cap_15_cis_window_200000_geno_filter_False_none_zscore_eqtl_input_" $surge_latent_factors_file $factor_pve_file $output_stem
fi

############################################
# Run SUSIE fine mapping
############################################
sample_names_file=$processed_pseudobulk_expression_dir"pseudobulk_scran_normalization_hvg_6000_regress_batch_True_individual_clustering_leiden_resolution_10.0_no_cap_15_sample_names.txt"
sh run_susie_finemapping_analysis.sh $surge_interaction_eqtl_dir $processed_genotype_dir $susie_input_data_dir $susie_results_dir $susie_visualization_dir $sample_names_file

############################################
# Check for overlap with coloc
############################################
if false; then
sh run_coloc_analysis.sh $surge_interaction_eqtl_dir $coloc_input_dir $coloc_dir $visualize_coloc_dir
fi


############################################
# Remap tfbs enrichment analysis
############################################
if false; then
sh run_remap_tfbs_enrichment_analysis.sh $surge_interaction_eqtl_dir $remap_tfbs_file $remap_tfbs_dir $latent_factor_interaction_eqtl_dir"latent_factor_interaction_hvg_6000_10.0_no_cap_15_cis_window_200000_geno_filter_False_none_zscore_eqtl_input_variant_gene_pairs.txt"
fi
############################################
# Check for overlap with gwas results
############################################
test_info_file=$latent_factor_interaction_eqtl_dir"latent_factor_interaction_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_input_variant_gene_pairs.txt"
if false; then
sh run_gwas_overlap_analysis.sh $surge_interaction_eqtl_dir $gwas_overlap_dir $coloc_dir"processed_gwas_studies.txt" $test_info_file $gwas_snp_info_dir
fi






############################################
# Run gene set enrichment analysis
############################################
model_loadings_file=$eqtl_factorization_results_dir"eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_factorization_vi_ard_heteroskedastic_results_k_init_10_seed_1_warmup_3000_ratio_variance_std_True_permute_False_lambda_1_round_geno_True_temper_U_S.txt"
model_pve_file=$eqtl_factorization_results_dir"eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_factorization_vi_ard_heteroskedastic_results_k_init_10_seed_1_warmup_3000_ratio_variance_std_True_permute_False_lambda_1_round_geno_True_temper_factor_pve.txt"
expr_file=$processed_pseudobulk_expression_dir"pseudobulk_scran_normalization_hvg_6000_regress_batch_True_individual_clustering_leiden_resolution_10.0_no_cap_15_none_sample_norm_zscore_gene_norm_normalized_expression.txt"
gene_names_file=$processed_pseudobulk_expression_dir"pseudobulk_scran_normalization_hvg_6000_regress_batch_True_individual_clustering_leiden_resolution_10.0_no_cap_15_gene_names.txt"
gene_set_enrichment_output_stem=$gene_set_enrichment_dir"standard_eqtl_hvg_6000_15_ard_heteroskedastic_rv_True_permute_False_seed_1_3000_warmup_lambda_1_v2_"
if false; then
sh gene_set_enrichment_analysis.sh $model_loadings_file $model_pve_file $expr_file $gene_names_file $gsea_data_dir $gene_set_enrichment_output_stem
fi

#############################################
## Visualize eqtl factorization
#############################################
if false; then
module load R/3.5.1
model_stem="eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_factorization_vi_ard_heteroskedastic_results_k_init_10_seed_1_warmup_3000_ratio_variance_std_True_permute_False_lambda_1_round_geno_True_temper_"
output_stem="standard_eqtl_hvg_6000_no_cap_15_ard_heteroskedastic_rv_True_permute_False_seed_1_3000_warmup_lambda_1_"
Rscript visualize_eqtl_factorization.R $processed_pseudobulk_expression_dir $eqtl_factorization_results_dir $eqtl_factorization_visualization_dir $model_stem $output_stem
fi






