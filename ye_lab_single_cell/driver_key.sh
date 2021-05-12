

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


######################
# Output directories
######################
# Root of all output directoreis
output_root="/work-zfs/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/"

# Directory containing processed genotype data
processed_genotype_dir=$output_root"processed_genotype/"

# Directory containing processed single cell expression
processed_expression_dir=$output_root"processed_expression/"

# Directory containing visualizations of processed single cell expression
visualize_processed_expression_dir=$output_root"visualize_processed_expression/"

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
sh preprocess_single_cell_expression.sh $input_h5py_file $processed_expression_dir $visualize_processed_expression_dir $gene_annotation_file $genotyped_individuals_file


######################
# Preprocess Genotype data
######################
# file containing list of individuals with have sc rna-seq for (generated by preprocess_single_cell_expression.sh)
sample_covariates_file=$processed_expression_dir"cluster_pseudobulk_leiden_no_cap_12_sample_covariates.txt"
if false; then
sh preprocess_genotype.sh $genotype_data_dir $processed_genotype_dir $sample_covariates_file $visualize_processed_genotype_dir
fi


######################
# Call latent-factor (PCA) interaction-eqtls
######################
if false; then
sh latent_factor_interaction_eqtl_driver_key.sh $processed_expression_dir $processed_genotype_dir $gene_annotation_file $latent_factor_interaction_eqtl_dir $visualize_latent_factor_interaction_eqtl_dir
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
sample_overlap_file=$eqtl_factorization_input_dir"eqtl_factorization_lf_interaction_eqtl_input_sample_overlap.txt"
expression_training_file=$eqtl_factorization_input_dir"eqtl_factorization_lf_interaction_eqtl_input_expression.npy"
genotype_training_file=$eqtl_factorization_input_dir"eqtl_factorization_lf_interaction_eqtl_input_genotype.npy"
covariate_file=$eqtl_factorization_input_dir"eqtl_factorization_lf_interaction_eqtl_input_covariates.txt"
num_latent_factors="20"
lambda_v="1"
model_name="eqtl_factorization_vi"
seed="1"

output_stem=$eqtl_factorization_results_dir"eqtl_factorization_results_lf_interaction_egenes_"$model_name"_results_k_init_"$num_latent_factors"_lambda_v_"$lambda_v"_seed_"$seed"_init4_"
if false; then
sbatch run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem
fi



#############################################
## Visualize eqtl factorization
#############################################
if false; then
module load R/3.5.1
Rscript visualize_eqtl_factorization.R $processed_expression_dir $eqtl_factorization_results_dir $eqtl_factorization_visualization_dir
fi


