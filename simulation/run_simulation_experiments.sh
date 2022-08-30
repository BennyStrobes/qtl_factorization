#!/bin/bash -l

#SBATCH
#SBATCH --time=50:00:00
#SBATCH --mem=20GB
#SBATCH --nodes=1




sim_data_dir="$1"
eqtl_results_dir="$2"



if false; then
missingness_fraction="0.1"
sbatch run_variance_explained_power_analysis.sh $sim_data_dir $eqtl_results_dir $missingness_fraction

missingness_fraction="0.3"
sbatch run_variance_explained_power_analysis.sh $sim_data_dir $eqtl_results_dir $missingness_fraction

missingness_fraction="0.5"
sbatch run_variance_explained_power_analysis.sh $sim_data_dir $eqtl_results_dir $missingness_fraction
fi

if false; then
sbatch run_num_component_selection_simulation_analysis.sh $sim_data_dir $eqtl_results_dir
fi