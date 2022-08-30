#!/bin/bash -l

#SBATCH
#SBATCH --time=50:00:00
#SBATCH --mem=20GB
#SBATCH --nodes=1





sim_data_dir="$1"
eqtl_results_dir="$2"
missingness_fraction="$3"







source ~/.bash_profile
conda activate surge

python run_variance_explained_power_analysis.py $sim_data_dir $eqtl_results_dir $missingness_fraction