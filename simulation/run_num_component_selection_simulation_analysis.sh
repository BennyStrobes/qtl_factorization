#!/bin/bash -l

#SBATCH
#SBATCH --time=50:00:00
#SBATCH --mem=20GB
#SBATCH --nodes=1





sim_data_dir="$1"
eqtl_results_dir="$2"
num_samples="$3"
t_statistic="$4"
missingness_fraction="$5"







source ~/.bash_profile
conda activate surge

python run_num_component_selection_simulation_analysis.py $sim_data_dir $eqtl_results_dir $num_samples $t_statistic $missingness_fraction