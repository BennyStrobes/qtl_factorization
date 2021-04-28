

######################
# Input data
######################
# Input single cell expression data (emailed by Meena Subramaniam on Nov. 18, 2019)
input_h5py_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/single_cell/input_data/CLUESImmVar_nonorm.V6.h5ad"

# Gene annotation file (hg19)
gene_annotation_file="/work-zfs/abattle4/lab_data/annotation/gencode.v19/gencode.v19.annotation.gene.txt"

# Directory containing genotype data
genotype_data_dir="/work-zfs/abattle4/lab_data/ye_lab_lupus_data_set/imputed_genotypes/"




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

# Directory containing standard eqtl dir
standard_eqtl_dir=$output_root"standard_eqtl/"

# Directory containing visualizations of standard eqtls
visualize_standard_eqtl_dir=$output_root"visualize_standard_eqtl/"

# Directory containing pre-processed eqtl factorization input files
eqtl_factorization_input_dir=$output_root"eqtl_factorization_input/"

# Directory containing eqtl factorization results
eqtl_factorization_results_dir=$output_root"eqtl_factorization_results/"

# Directory containing visualizations of eqtl factorization results
eqtl_factorization_visualization_dir=$output_root"visualize_eqtl_factorization/"







######################
# Preprocess Genotype data
######################
if false; then
sbatch preprocess_genotype.sh $genotype_data_dir $processed_genotype_dir
fi


######################
# Preprocess single cell expression
######################
if false; then
sbatch preprocess_single_cell_expression.sh $input_h5py_file $processed_expression_dir $visualize_processed_expression_dir $gene_annotation_file
fi