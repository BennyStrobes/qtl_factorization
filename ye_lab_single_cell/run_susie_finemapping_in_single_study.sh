#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


sig_eqtl_file="$1"
all_eqtl_file="$2"
eqtl_study_name="$3"
processed_genotype_dir="$4"
output_processed_data_root="$5"
output_results_root="$6"
output_visualization_root="$7"
surge_latent_factor_file="$8"
surge_latent_factor_num="$9"
sample_names_file="${10}"


echo "preparing susie data"
if false; then
python prepare_data_for_susie_finemapping.py $sig_eqtl_file $all_eqtl_file $processed_genotype_dir $output_processed_data_root $surge_latent_factor_file $surge_latent_factor_num $sample_names_file
fi


if false; then
module load R/3.5.1
Rscript run_susie_finemapping.R $output_processed_data_root $output_results_root $output_visualization_root
fi