#!/bin/bash -l

#SBATCH
#SBATCH --time=40:00:00
#SBATCH --partition=shared
#SBATCH --mem=40GB

#SBATCH --nodes=1




expression_training_file="$1"
genotype_training_file="$2"
covariate_file="$3"
sample_overlap_file="$4"
num_latent_factors="$5"
lambda_v="$6"
model_name="$7"
seed="$8"
output_stem="$9"
variance_param="${10}"

module load python/3.7-anaconda

python run_eqtl_factorization.py $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param