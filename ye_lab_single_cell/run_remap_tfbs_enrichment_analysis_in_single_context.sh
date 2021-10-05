#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --mem=50GB
#SBATCH --nodes=1



surge_interaction_eqtl_dir="$1"
remap_tfbs_file="$2"
remap_tfbs_dir="$3"
variant_info_file="$4"
gsea_data_dir="$5"
latent_context="$6"


echo $latent_context

surge_interaction_sig_eqtl_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_interaction_eqtl_results_latent_factor_"$latent_context"_genome_wide_signficant_bf_fdr_0.05.txt"
output_file=$remap_tfbs_dir"surge_interaction_"$latent_context"_remap_overlap.txt"
gsea_test_genes_file=$remap_tfbs_dir"surge_interaction_"$latent_context"_remap_test_genes.txt"
gsea_bgrd_genes_file=$remap_tfbs_dir"surge_interaction_"$latent_context"_remap_bgrd_genes.txt"
gsea_output_root=$remap_tfbs_dir"surge_interaction_"$latent_context"_remap_enriched_tf_gsea_"
python run_remap_tfbs_enrichment_analysis.py $surge_interaction_sig_eqtl_file $remap_tfbs_file $output_file $variant_info_file $gsea_test_genes_file $gsea_bgrd_genes_file

source ~/.bash_profile

module load python/2.7-anaconda53
gsea $gsea_test_genes_file $gsea_bgrd_genes_file $gsea_data_dir"h.all.v5.1.symbols.gmt.txt" $gsea_output_root"enrichments_h_all_v5.txt"
python order_and_threshold_gsea.py $gsea_output_root"enrichments_h_all_v5.txt"	
gsea $gsea_test_genes_file $gsea_bgrd_genes_file $gsea_data_dir"c2.cp.kegg.v5.1.symbols.gmt.txt" $gsea_output_root"enrichments_c2_cp_kegg.txt"
python order_and_threshold_gsea.py $gsea_output_root"enrichments_c2_cp_kegg.txt"
gsea $gsea_test_genes_file $gsea_bgrd_genes_file $gsea_data_dir"c5.bp.v5.1.symbols.gmt.txt" $gsea_output_root"enrichments_c5_bp.txt"
python order_and_threshold_gsea.py $gsea_output_root"enrichments_c5_bp.txt"
gsea $gsea_test_genes_file $gsea_bgrd_genes_file $gsea_data_dir"c5.mf.v5.1.symbols.gmt.txt" $gsea_output_root"enrichments_c5_mf.txt"
python order_and_threshold_gsea.py $gsea_output_root"enrichments_c5_mf.txt"
gsea $gsea_test_genes_file $gsea_bgrd_genes_file $gsea_data_dir"c2.cp.biocarta.v5.1.symbols.gmt.txt" $gsea_output_root"enrichments_c2_cp_biocarta.txt"
python order_and_threshold_gsea.py $gsea_output_root"enrichments_c2_cp_biocarta.txt"	
gsea $gsea_test_genes_file $gsea_bgrd_genes_file $gsea_data_dir"c2.cp.reactome.v5.1.symbols.gmt.txt" $gsea_output_root"enrichments_c2_cp_reactome.txt"
python order_and_threshold_gsea.py $gsea_output_root"enrichments_c2_cp_reactome.txt"	
gsea $gsea_test_genes_file $gsea_bgrd_genes_file $gsea_data_dir"c4.cgn.v5.1.symbols.gmt.txt" $gsea_output_root"enrichments_c4_cgn.txt"
python order_and_threshold_gsea.py $gsea_output_root"enrichments_c4_cgn.txt"
gsea $gsea_test_genes_file $gsea_bgrd_genes_file $gsea_data_dir"c4.cm.v5.1.symbols.gmt.txt" $gsea_output_root"enrichments_c4_cm.txt"
python order_and_threshold_gsea.py $gsea_output_root"enrichments_c4_cm.txt"
gsea $gsea_test_genes_file $gsea_bgrd_genes_file $gsea_data_dir"c6.all.v5.1.symbols.gmt.txt" $gsea_output_root"enrichments_c6_all.txt"
python order_and_threshold_gsea.py $gsea_output_root"enrichments_c6_all.txt"
gsea $gsea_test_genes_file $gsea_bgrd_genes_file $gsea_data_dir"c7.all.v5.1.symbols.gmt.txt" $gsea_output_root"enrichments_c7_all.txt"
python order_and_threshold_gsea.py $gsea_output_root"enrichments_c7_all.txt"
gsea $gsea_test_genes_file $gsea_bgrd_genes_file $gsea_data_dir"human_genesigdb.gmt" $gsea_output_root"enrichments_human_genesigdb.txt"
python order_and_threshold_gsea.py $gsea_output_root"enrichments_human_genesigdb.txt"





