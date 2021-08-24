################################
# Input data
#################################
# File containing raw, log2(CPM+1) data
normalized_expression_file="/work-zfs/abattle4/lab_data/sc_endo_diff/counts.tsv"

# File containing meta-data for each cell
meta_data_file="/work-zfs/abattle4/lab_data/sc_endo_diff/cell_metadata_cols.tsv"

# File containing vcf files for each individual
genotype_dir="/work-zfs/abattle4/prashanthi/sc-endo/data/genotypes/"

# Gencode hg19 gene annotation file
gene_annotation_file="/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt"

# Genotype PC file
genotype_pc_file="/work-zfs/abattle4/prashanthi/sc-endo/data/genotypes/genotype_PCs.txt"

# File generated from Cuomo et al group which contains info describing cell state of each cell
cell_state_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell_differentiation_cuomo_data_ase/input_data/sce_merged_afterqc_filt_allexpts_pseudotimeandmodules_20180618.tsv"


################################
# Output Directories
################################
# Root directory
root_directory="/work-zfs/abattle4/bstrober/qtl_factorization/cuomo_differentiation_single_cell/"
# Directory containing pre-processed genotype 
processed_genotype_dir=$root_directory"processed_genotype/"
# Directory containing pre-processed expression
processed_expression_dir=$root_directory"processed_expression/"
# Directory containing visualizations of processed gene expression
visualize_processed_expression_dir=$root_directory"visualize_processed_expression/"
# Directory containing standard eqtl input data
standard_eqtl_input_dir=$root_directory"standard_eqtl_input_data/"
# Directory containing standard eqtl results
standard_eqtl_results_dir=$root_directory"standard_eqtl_results/"
# Directory containing eQTL factorization input
eqtl_factorization_input_dir=$root_directory"eqtl_factorization_input/"
# Directory containing eqtl factorization results
eqtl_factorization_results_dir=$root_directory"eqtl_factorization_results/"
# Directory containing visualizations of eqtl factorization results
eqtl_factorization_visualization_dir=$root_directory"visualize_eqtl_factorization/"


################################
# Preprocess genotype data
#################################
if false; then
sh preprocess_genotype.sh $genotype_dir $processed_genotype_dir
fi


################################
# Preprocess gene expression data
#################################
if false; then
sbatch preprocess_expression.sh $normalized_expression_file $meta_data_file $processed_genotype_dir $gene_annotation_file $genotype_pc_file $cell_state_file $processed_expression_dir $visualize_processed_expression_dir
fi



################################
# Run standard eQTL analysis
#################################
if false; then
sh run_standard_eqtl_analysis.sh $processed_expression_dir $processed_genotype_dir $gene_annotation_file $standard_eqtl_input_dir $standard_eqtl_results_dir
fi

################################
# Prepare input data for eQTL factorization
#################################
if false; then
sh preprocess_data_for_eqtl_factorization.sh $standard_eqtl_input_dir $standard_eqtl_results_dir $eqtl_factorization_input_dir $processed_expression_dir
fi

#############################################
## Run eqtl factorization!
#############################################
input_data_stem="no_repeats_no_hla_eqtl_factorization_standard_eqtl_scanpy_4000_hvg"
sample_overlap_file=$eqtl_factorization_input_dir$input_data_stem"_sample_overlap.txt"
expression_training_file=$eqtl_factorization_input_dir$input_data_stem"_expression.npy"
genotype_training_file=$eqtl_factorization_input_dir$input_data_stem"_unnormalized_genotype.npy"
covariate_file=$eqtl_factorization_input_dir$input_data_stem"_covariates.txt"
num_latent_factors="5"
lambda_v="1"
variance_param="1e-3"
ard_variance_param="1e-16"
seed="2"
model_name="eqtl_factorization_vi_ard"
ratio_variance_standardization="True"
permutation_type="False"
warmup_iterations="5"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_remove_heterozygous_outliers_"
if false; then
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations

seed="2"
model_name="eqtl_factorization_vi_ard_full_component_update"
ratio_variance_standardization="True"
permutation_type="False"
warmup_iterations="5"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_remove_heterozygous_outliers_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations
fi

if false; then
seed="2"
model_name="eqtl_factorization_vi_ard_heteroskedastic"
ratio_variance_standardization="True"
permutation_type="False"
warmup_iterations="5"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations

seed="2"
model_name="eqtl_factorization_vi_no_factorization"
ratio_variance_standardization="True"
permutation_type="False"
warmup_iterations="5"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations

seed="2"
model_name="eqtl_factorization_vi_ard_full_component_update"
ratio_variance_standardization="True"
permutation_type="False"
warmup_iterations="5"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations


seed="2"
model_name="eqtl_factorization_vi_ard_heteroskedastic_full_component_update"
ratio_variance_standardization="True"
permutation_type="False"
warmup_iterations="5"
output_stem=$eqtl_factorization_results_dir$input_data_stem"_"$model_name"_results_k_init_"$num_latent_factors"_seed_"$seed"_warmup_"$warmup_iterations"_ratio_variance_std_"$ratio_variance_standardization"_permute_"$permutation_type"_"
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations
fi


if false; then
module load R/3.5.1
model_stem="no_repeats_no_hla_eqtl_factorization_standard_eqtl_scanpy_4000_hvg_eqtl_factorization_vi_ard_results_k_init_5_seed_2_warmup_5_ratio_variance_std_True_permute_False_remove_heterozygous_outliers_temper_"
output_stem="no_repeats_no_hla_standard_eqtl_remove_heterozygous_outliers_scanpy_4000_hvg_ard_rv_True_permute_False_seed_1_k_5_warmup_5_"
Rscript visualize_eqtl_factorization.R $processed_expression_dir $eqtl_factorization_results_dir $eqtl_factorization_visualization_dir $model_stem $output_stem

module load R/3.5.1
model_stem="no_repeats_no_hla_eqtl_factorization_standard_eqtl_scanpy_4000_hvg_eqtl_factorization_vi_ard_full_component_update_results_k_init_5_seed_2_warmup_5_ratio_variance_std_True_permute_False_remove_heterozygous_outliers_temper_"
output_stem="no_repeats_no_hla_standard_eqtl_full_component_update_remove_heterozygous_outliers_scanpy_4000_hvg_ard_rv_True_permute_False_seed_1_k_5_warmup_5_"
Rscript visualize_eqtl_factorization.R $processed_expression_dir $eqtl_factorization_results_dir $eqtl_factorization_visualization_dir $model_stem $output_stem

module load R/3.5.1
model_stem="no_repeats_no_hla_eqtl_factorization_standard_eqtl_scanpy_4000_hvg_eqtl_factorization_vi_ard_results_k_init_5_seed_2_warmup_5_ratio_variance_std_True_permute_False_temper_"
output_stem="no_repeats_no_hla_standard_eqtl_scanpy_4000_hvg_ard_rv_True_permute_False_seed_1_k_5_warmup_5_"
Rscript visualize_eqtl_factorization.R $processed_expression_dir $eqtl_factorization_results_dir $eqtl_factorization_visualization_dir $model_stem $output_stem
fi

module load R/3.5.1
model_stem="no_repeats_no_hla_eqtl_factorization_standard_eqtl_scanpy_4000_hvg_eqtl_factorization_vi_ard_full_component_update_results_k_init_5_seed_2_warmup_5_ratio_variance_std_True_permute_False_temper_"
output_stem="no_repeats_no_hla_standard_eqtl_full_component_update_scanpy_4000_hvg_ard_rv_True_permute_False_seed_1_k_5_warmup_5_"
Rscript visualize_eqtl_factorization.R $processed_expression_dir $eqtl_factorization_results_dir $eqtl_factorization_visualization_dir $model_stem $output_stem

if false; then
module load R/3.5.1
model_stem="no_repeats_no_hla_eqtl_factorization_standard_eqtl_scanpy_4000_hvg_eqtl_factorization_vi_ard_heteroskedastic_full_component_update_results_k_init_5_seed_2_warmup_5_ratio_variance_std_True_permute_False_temper_"
output_stem="no_repeats_no_hla_standard_eqtl_heteroskedastic_full_component_update_scanpy_4000_hvg_ard_rv_True_permute_False_seed_1_k_5_warmup_5_"
Rscript visualize_eqtl_factorization.R $processed_expression_dir $eqtl_factorization_results_dir $eqtl_factorization_visualization_dir $model_stem $output_stem

module load R/3.5.1
model_stem="no_repeats_no_hla_eqtl_factorization_standard_eqtl_scanpy_4000_hvg_eqtl_factorization_vi_ard_heteroskedastic_results_k_init_5_seed_2_warmup_5_ratio_variance_std_True_permute_False_temper_"
output_stem="no_repeats_no_hla_standard_eqtl_heteroskedastic_full_scanpy_4000_hvg_ard_rv_True_permute_False_seed_1_k_5_warmup_5_"
Rscript visualize_eqtl_factorization.R $processed_expression_dir $eqtl_factorization_results_dir $eqtl_factorization_visualization_dir $model_stem $output_stem
fi



