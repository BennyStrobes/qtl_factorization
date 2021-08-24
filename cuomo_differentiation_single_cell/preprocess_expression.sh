#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=lrgmem
#SBATCH --mem=80GB

#SBATCH --nodes=1


normalized_expression_file="$1"
meta_data_file="$2"
processed_genotype_dir="$3"
gene_annotation_file="$4"
genotype_pc_file="$5"
cell_state_file="$6"
processed_expression_dir="$7"
visualize_processed_expression_dir="$8"

module load R/3.6.1
module load python/3.7.4-anaconda
python preprocess_expression.py $normalized_expression_file $meta_data_file $processed_genotype_dir $gene_annotation_file $genotype_pc_file $cell_state_file $processed_expression_dir


if false; then
module load R/3.5.1
Rscript visualize_processed_expression.R $processed_expression_dir $visualize_processed_expression_dir
fi