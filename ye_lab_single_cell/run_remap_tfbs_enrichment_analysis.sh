#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


surge_interaction_eqtl_dir="$1"
remap_tfbs_file="$2"
remap_tfbs_dir="$3"
variant_info_file="$4"


latent_context="4"
surge_interaction_sig_eqtl_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_interaction_eqtl_results_latent_factor_"$latent_context"_genome_wide_signficant_bf_fdr_0.05.txt"
output_file=$remap_tfbs_dir"surge_interaction_"$latent_context"_remap_overlap.txt"

python run_remap_tfbs_enrichment_analysis.py $surge_interaction_sig_eqtl_file $remap_tfbs_file $output_file $variant_info_file