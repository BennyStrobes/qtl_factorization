#!/bin/bash -l

#SBATCH
#SBATCH --time=30:00:00
#SBATCH --mem=600GB
#SBATCH --partition=lrgmem
#SBATCH --nodes=1

input_h5py_file="$1"
processed_expression_dir="$2"
visualize_processed_expression_dir="$3"
gene_annotation_file="$4"
genotyped_individuals_file="$5"
processed_pseudobulk_expression_dir="$6"
visualize_processed_pseudobulk_expression_dir="$7"
regress_out_batch="$8"
isg_score_file="$9"
cell_isg_score_file="${10}"

module load R/3.6.1
module load python/3.7.4-anaconda


# Preprocess single cell expression
if false; then
python preprocess_scran_single_cell_expression.py $input_h5py_file $processed_expression_dir $gene_annotation_file $genotyped_individuals_file $regress_out_batch
fi

# Visualize preprocessed single cell expression
module load R/3.5.1
if false; then
Rscript visualize_scran_processed_single_cell_expression.R $processed_expression_dir $visualize_processed_expression_dir $regress_out_batch
fi

cluster_resolution="10.0"
cluster_type="no_cap"
min_cells_per_cluster="15"
num_hvg="6000"
module load R/3.6.1
module load python/3.7.4-anaconda
# Generate pseudobulk expression
python generate_pseudobulk_expression.py $processed_expression_dir $processed_pseudobulk_expression_dir $genotyped_individuals_file $cluster_resolution $cluster_type $regress_out_batch $num_hvg $isg_score_file $cell_isg_score_file $min_cells_per_cluster

if false; then
cluster_resolution="10.0"
num_hvg="4000"
module load R/3.6.1
module load python/3.7.4-anaconda
# Generate pseudobulk expression
python generate_pseudobulk_expression.py $processed_expression_dir $processed_pseudobulk_expression_dir $genotyped_individuals_file $cluster_resolution $cluster_type $regress_out_batch $num_hvg $isg_score_file $cell_isg_score_file

cluster_resolution="10.0"
num_hvg="2000"
module load R/3.6.1
module load python/3.7.4-anaconda
# Generate pseudobulk expression
python generate_pseudobulk_expression.py $processed_expression_dir $processed_pseudobulk_expression_dir $genotyped_individuals_file $cluster_resolution $cluster_type $regress_out_batch $num_hvg $isg_score_file $cell_isg_score_file
fi


if false; then
module load R/3.5.1
sample_level_normalization="qn"
gene_level_normalization="zscore"
Rscript visualize_pseudobulk_expression.R $processed_pseudobulk_expression_dir $cluster_resolution $visualize_processed_pseudobulk_expression_dir $regress_out_batch $gene_level_normalization $sample_level_normalization

sample_level_normalization="qn"
gene_level_normalization="ign"
Rscript visualize_pseudobulk_expression.R $processed_pseudobulk_expression_dir $cluster_resolution $visualize_processed_pseudobulk_expression_dir $regress_out_batch $gene_level_normalization $sample_level_normalization

sample_level_normalization="none"
gene_level_normalization="zscore"
Rscript visualize_pseudobulk_expression.R $processed_pseudobulk_expression_dir $cluster_resolution $visualize_processed_pseudobulk_expression_dir $regress_out_batch $gene_level_normalization $sample_level_normalization

sample_level_normalization="none"
gene_level_normalization="ign"
Rscript visualize_pseudobulk_expression.R $processed_pseudobulk_expression_dir $cluster_resolution $visualize_processed_pseudobulk_expression_dir $regress_out_batch $gene_level_normalization $sample_level_normalization

fi
























if false; then
## CURRENTLY OLD
python preprocess_single_cell_expression.py $input_h5py_file $processed_expression_dir $gene_annotation_file $genotyped_individuals_file
fi

if false; then
module load R/3.5.1
Rscript visualize_processed_single_cell_expression.R $processed_expression_dir $visualize_processed_expression_dir
fi
