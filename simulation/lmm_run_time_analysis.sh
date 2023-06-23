#!/bin/bash -l

#SBATCH
#SBATCH --time=50:00:00
#SBATCH --mem=35GB
#SBATCH --nodes=1





sim_data_dir="$1"
eqtl_results_dir="$2"
run_time_iter="$3"


module load r/3.6.3

Rscript run_lmm_run_time_analysis.R $sim_data_dir $eqtl_results_dir $run_time_iter