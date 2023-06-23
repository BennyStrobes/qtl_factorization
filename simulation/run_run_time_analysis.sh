#!/bin/bash -l

#SBATCH
#SBATCH --time=72:00:00
#SBATCH --mem=35GB
#SBATCH --nodes=1





sim_data_dir="$1"
eqtl_results_dir="$2"
run_time_iter="$3"





source ~/.bash_profile
conda activate surge

python run_run_time_analysis.py $sim_data_dir $eqtl_results_dir $run_time_iter
