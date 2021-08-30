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
sample_covariates_file=$processed_pseudobulk_expression_dir"pseudobulk_scran_normalization_hvg_6000_regress_batch_True_individual_clustering_leiden_resolution_10.0_no_cap_sample_covariates.txt"
if false; then
sh preprocess_genotype.sh $genotype_data_dir $processed_genotype_dir $sample_covariates_file $visualize_processed_genotype_dir
fi


######################
# Call latent-factor (PCA) interaction-eqtls
######################
sh latent_factor_interaction_eqtl_driver_key.sh $processed_expression_dir $processed_pseudobulk_expression_dir $processed_genotype_dir $gene_annotation_file $latent_factor_interaction_eqtl_dir $visualize_latent_factor_interaction_eqtl_dir


#############################################
# Prepare data for eQTL factorization 
#############################################
if false; then
sh preprocess_data_for_eqtl_factorization.sh $latent_factor_interaction_eqtl_dir $eqtl_factorization_input_dir
fi



#############################################
## Run eqtl factorization!
#############################################
input_data_stem="eqtl_factorization_standard_eqtl_10.0_none_zscore_capped"
sample_overlap_file=$eqtl_factorization_input_dir$input_data_stem"_eqtl_input_sample_overlap.txt"
expression_training_file=$eqtl_factorization_input_dir$input_data_stem"_eqtl_input_expression.npy"
genotype_training_file=$eqtl_factorization_input_dir$input_data_stem"_eqtl_input_unnormalized_genotype.npy"
covariate_file=$eqtl_factorization_input_dir$input_data_stem"_eqtl_input_covariates.txt"
num_latent_factors="10"
lambda_v="1"
variance_param="1e-3"
ard_variance_param="1e-16"
seed="1"
permutation_type="False"
model_name="eqtl_factorization_vi_ard_heteroskedastic"
ratio_variance_standardization="False"
warmup_iterations="5"

output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_lambda_"$lambda_v"_round_genotype_"
if false; then
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations
fi



model_name="eqtl_factorization_vi_no_factorization"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_"
if false; then
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations
fi


if false; then

seed="1"
model_name="eqtl_factorization_vi_ard"
ratio_variance_standardization="True"
permutation_type="False"
warmup_iterations="6"


output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations


seed="1"
model_name="eqtl_factorization_vi_ard"
ratio_variance_standardization="False"
permutation_type="False"
warmup_iterations="0"


output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations

seed="1"
model_name="eqtl_factorization_vi_ard"
ratio_variance_standardization="False"
permutation_type="False"
warmup_iterations="6"


output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations


seed="1"
model_name="eqtl_factorization_vi_ard"
ratio_variance_standardization="True"
permutation_type="interaction_only"
warmup_iterations="0"


output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations

seed="1"
model_name="eqtl_factorization_vi_ard"
ratio_variance_standardization="True"
permutation_type="fixed_and_interaction"
warmup_iterations="0"


output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations
fi

#############################################
## Visualize eqtl factorization
#############################################
if false; then
module load R/3.5.1
model_stem="eqtl_factorization_standard_eqtl_10.0_none_zscore_capped_eqtl_factorization_vi_ard_heteroskedastic_results_k_init_10_seed_1_warmup_0_ratio_variance_std_True_permute_False_round_genotype_temper_"
output_stem="standard_eqtl_heteroskedastic_rv_True_permute_False_seed_1_0_warmup_lambda_1_"
Rscript visualize_eqtl_factorization.R $processed_pseudobulk_expression_dir $eqtl_factorization_results_dir $eqtl_factorization_visualization_dir $model_stem $output_stem
fi




if false; then
module load R/3.5.1
Rscript visualize_eqtl_lda.R $processed_pseudobulk_expression_dir $eqtl_factorization_results_dir $eqtl_factorization_visualization_dir
fi