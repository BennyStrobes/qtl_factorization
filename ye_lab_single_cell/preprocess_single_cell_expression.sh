#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --mem=180GB
#SBATCH --partition=lrgmem
#SBATCH --nodes=1

input_h5py_file="$1"
processed_expression_dir="$2"
visualize_processed_expression_dir="$3"
gene_annotation_file="$4"
genotype_id_mapping_file="$5"


module load python/3.7.4-anaconda

min_fraction_cells=".05"
min_genes="400"
transformation_type="log_transform"
python preprocess_single_cell_expression.py $input_h5py_file $processed_expression_dir $gene_annotation_file $min_fraction_cells $min_genes $transformation_type $genotype_id_mapping_file


if false; then
module load R/3.5.1
Rscript visualize_processed_single_cell_expression.R $processed_expression_dir $visualize_processed_expression_dir
fi



