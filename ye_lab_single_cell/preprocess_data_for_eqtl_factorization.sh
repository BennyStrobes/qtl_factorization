#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --mem=200GB
#SBATCH --partition=bigmem
#SBATCH -A abattle44_bigmem
#SBATCH --nodes=1



latent_factor_interaction_eqtl_dir="$1"
eqtl_factorization_input_dir="$2"




num_genes="2000"

echo "STARTING"
python preprocess_data_for_eqtl_factorization_based_on_standard_eqtls.py $num_genes $latent_factor_interaction_eqtl_dir $eqtl_factorization_input_dir



if false; then
python preprocess_data_for_eqtl_factorization_based_on_random_variant_per_gene.py $num_genes $latent_factor_interaction_eqtl_dir $eqtl_factorization_input_dir
fi