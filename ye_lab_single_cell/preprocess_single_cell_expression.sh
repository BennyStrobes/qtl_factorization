#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --mem=400GB
#SBATCH --partition=lrgmem
#SBATCH --nodes=1

input_h5py_file="$1"
processed_expression_dir="$2"
visualize_processed_expression_dir="$3"
gene_annotation_file="$4"
genotyped_individuals_file="$5"

module load R/3.5.1
module load python/3.7.4-anaconda
python preprocess_single_cell_expression.py $input_h5py_file $processed_expression_dir $gene_annotation_file $genotyped_individuals_file


if false; then
module load R/3.5.1
Rscript visualize_processed_single_cell_expression.R $processed_expression_dir $visualize_processed_expression_dir
fi

