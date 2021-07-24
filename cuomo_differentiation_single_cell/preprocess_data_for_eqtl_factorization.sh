
#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --mem=250GB
#SBATCH --partition=lrgmem
#SBATCH --nodes=1


standard_eqtl_input_dir="$1"
standard_eqtl_results_dir="$2"
eqtl_factorization_input_dir="$3"


python preprocess_data_for_eqtl_factorization_based_on_standard_eqtls.py $standard_eqtl_input_dir $standard_eqtl_results_dir $eqtl_factorization_input_dir
