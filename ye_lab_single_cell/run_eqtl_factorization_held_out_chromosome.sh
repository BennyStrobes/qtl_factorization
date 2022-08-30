#!/bin/bash -l

#SBATCH
#SBATCH --time=72:00:00
#SBATCH --mem=50GB
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
ard_variance_param="${11}"
ratio_variance_standardization="${12}"
permutation_type="${13}"
warmup_iterations="${14}"
round_genotype="${15}"
data_filter="${16}"
test_names_file="${17}"
delta_elbo_thresh="${18}"
held_out_chrom="${19}"

source ~/.bash_profile
conda activate surge 

python run_eqtl_factorization_held_out_chromosome.py $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem $variance_param $ard_variance_param $ratio_variance_standardization $permutation_type $warmup_iterations $round_genotype $data_filter $test_names_file $delta_elbo_thresh $held_out_chrom