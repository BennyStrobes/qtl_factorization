#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --mem=250GB
#SBATCH --partition=lrgmem
#SBATCH --nodes=1



latent_factor_interaction_eqtl_dir="$1"
eqtl_factorization_input_dir="$2"




num_genes="900"




python preprocess_data_for_eqtl_factorization_based_on_latent_factor_interaction_eqtls.py $num_genes $latent_factor_interaction_eqtl_dir $eqtl_factorization_input_dir

if false; then
python preprocess_data_for_eqtl_factorization_based_on_random_variant_per_gene.py $num_genes $latent_factor_interaction_eqtl_dir $eqtl_factorization_input_dir
fi