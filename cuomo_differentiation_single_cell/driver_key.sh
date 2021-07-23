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
sh preprocess_expression.sh $normalized_expression_file $meta_data_file $processed_genotype_dir $gene_annotation_file $genotype_pc_file $cell_state_file $processed_expression_dir $visualize_processed_expression_dir
fi




