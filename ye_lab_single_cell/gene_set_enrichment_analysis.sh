#!/bin/bash -l

#SBATCH
#SBATCH --time=50:00:00
#SBATCH --partition=shared
#SBATCH --mem=40GB

#SBATCH --nodes=1


model_loadings_file="$1"
expr_file="$2"
gene_names_file="$3"
gsea_data_dir="$4"
gene_set_enrichment_output_stem="$5"

num_components="10"

for component_num in $(seq 0 $(($num_components-1))); do 
	component_num="5"
	echo $component_num
	component_output_stem=$gene_set_enrichment_output_stem"component_"$component_num"_"
	python prepare_gsea_input_data.py $model_loadings_file $expr_file $gene_names_file $component_num $component_output_stem


	source ~/.bash_profile

	module load python/2.7-anaconda53
	gsea $component_output_stem"test_genes.txt" $component_output_stem"background_genes.txt" $gsea_data_dir"h.all.v5.1.symbols.gmt.txt" $component_output_stem"enrichments_h_all_v5.txt"
done



