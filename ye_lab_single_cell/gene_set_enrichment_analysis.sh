#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --mem=10GB

#SBATCH --nodes=1


model_loadings_file="$1"
model_pve_file="$2"
expr_file="$3"
gene_names_file="$4"
gsea_data_dir="$5"
gene_set_enrichment_output_stem="$6"



num_components="10"
for component_num in $(seq 3 $(($num_components-1))); do 
	component_output_stem=$gene_set_enrichment_output_stem"component_"$component_num"_"
	python prepare_gsea_input_data.py $model_loadings_file $model_pve_file $expr_file $gene_names_file $component_num $component_output_stem

	source ~/.bash_profile

	module load python/2.7-anaconda53
	gsea $component_output_stem"test_genes.txt" $component_output_stem"background_genes.txt" $gsea_data_dir"h.all.v5.1.symbols.gmt.txt" $component_output_stem"enrichments_h_all_v5.txt"
	python order_and_threshold_gsea.py $component_output_stem"enrichments_h_all_v5.txt"
	gsea $component_output_stem"test_genes.txt" $component_output_stem"background_genes.txt" $gsea_data_dir"c2.cp.kegg.v5.1.symbols.gmt.txt" $component_output_stem"enrichments_c2_cp_kegg.txt"
	python order_and_threshold_gsea.py $component_output_stem"enrichments_c2_cp_kegg.txt"
	gsea $component_output_stem"test_genes.txt" $component_output_stem"background_genes.txt" $gsea_data_dir"c5.bp.v5.1.symbols.gmt.txt" $component_output_stem"enrichments_c5_bp.txt"
	python order_and_threshold_gsea.py $component_output_stem"enrichments_c5_bp.txt"
done


if false; then
python organize_gsea_results_cross_component.py $gene_set_enrichment_output_stem "enrichments_c5_bp.txt"

python organize_gsea_results_cross_component.py $gene_set_enrichment_output_stem "enrichments_c2_cp_kegg.txt"
fi
if false; then

module load R/3.5.1
Rscript visualize_gsea_results.R $gene_set_enrichment_output_stem
fi