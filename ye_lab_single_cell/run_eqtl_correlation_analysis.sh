#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1




surge_interaction_eqtl_dir="$1"
eqtl_correlation_dir="$2"
python extract_eqtl_correlations.py $surge_interaction_eqtl_dir $eqtl_correlation_dir


module load R/3.5.1
Rscript visualize_eqtl_correlations.R $eqtl_correlation_dir
