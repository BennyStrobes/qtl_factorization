#!/bin/bash -l

#SBATCH
#SBATCH --time=40:00:00
#SBATCH --partition=shared
#SBATCH --mem=50GB
#SBATCH --nodes=1


surge_interaction_eqtl_dir="$1"
remap_tfbs_file="$2"
remap_tfbs_dir="$3"
variant_info_file="$4"
gsea_data_dir="$5"

for latent_context in $(seq 1 $((10))); do 
	sh run_remap_tfbs_enrichment_analysis_in_single_context.sh $surge_interaction_eqtl_dir $remap_tfbs_file $remap_tfbs_dir $variant_info_file $gsea_data_dir $latent_context 
done





