#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1



surge_interaction_eqtl_dir="$1"
coloc_input_dir="$2"
coloc_dir="$3"
visualize_coloc_dir="$4"



processed_gwas_studies_file=$coloc_dir"processed_gwas_studies.txt"
if false; then
python prepare_gwas_data_for_coloc_analysis.py $coloc_input_dir $processed_gwas_studies_file $coloc_dir
fi


eqtl_study_name="surge_latent_factor_1_interaction"
sig_eqtl_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_factor_1_interaction_eqtl_results_genome_wide_signficant_bf_fdr_0.05.txt"
all_eqtl_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_factor_1_interaction_eqtl_results_merged.txt"
if false; then
while read gwas_study_name gwas_study_file_root; do
	echo $gwas_study_name
	output_root=$coloc_dir$eqtl_study_name"_"$gwas_study_name"_coloc_"
	sh run_coloc_for_study_pairing.sh $sig_eqtl_file $all_eqtl_file $gwas_study_file_root $output_root $eqtl_study_name $gwas_study_name "version1"
done<$processed_gwas_studies_file
fi




eqtl_study_name="surge_latent_factor_2_interaction"
sig_eqtl_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_factor_2_interaction_eqtl_results_genome_wide_signficant_bf_fdr_0.05.txt"
all_eqtl_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_factor_2_interaction_eqtl_results_merged.txt"
if false; then
while read gwas_study_name gwas_study_file_root; do
	echo $gwas_study_name
	output_root=$coloc_dir$eqtl_study_name"_"$gwas_study_name"_coloc_"
	sh run_coloc_for_study_pairing.sh $sig_eqtl_file $all_eqtl_file $gwas_study_file_root $output_root $eqtl_study_name $gwas_study_name "version1"
done<$processed_gwas_studies_file
fi

eqtl_study_name="surge_latent_factor_3_interaction"
sig_eqtl_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_factor_3_interaction_eqtl_results_genome_wide_signficant_bf_fdr_0.05.txt"
all_eqtl_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_factor_3_interaction_eqtl_results_merged.txt"
if false; then
while read gwas_study_name gwas_study_file_root; do
	echo $gwas_study_name
	output_root=$coloc_dir$eqtl_study_name"_"$gwas_study_name"_coloc_"
	sh run_coloc_for_study_pairing.sh $sig_eqtl_file $all_eqtl_file $gwas_study_file_root $output_root $eqtl_study_name $gwas_study_name "version1"
done<$processed_gwas_studies_file
fi



eqtl_study_name="surge_latent_factor_4_interaction"
sig_eqtl_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_factor_4_interaction_eqtl_results_genome_wide_signficant_bf_fdr_0.05.txt"
all_eqtl_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_factor_4_interaction_eqtl_results_merged.txt"
if false; then
while read gwas_study_name gwas_study_file_root; do
	echo $gwas_study_name
	output_root=$coloc_dir$eqtl_study_name"_"$gwas_study_name"_coloc_"
	sh run_coloc_for_study_pairing.sh $sig_eqtl_file $all_eqtl_file $gwas_study_file_root $output_root $eqtl_study_name $gwas_study_name "version1"
done<$processed_gwas_studies_file
fi

eqtl_study_name="surge_latent_factor_4_interaction_big_window"
sig_eqtl_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_4_interaction_eqtl_results_genome_wide_signficant_bf_fdr_0.05.txt"
all_eqtl_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_4_interaction_eqtl_results_merged.txt"
if false; then
while read gwas_study_name gwas_study_file_root; do
	echo $gwas_study_name
	output_root=$coloc_dir$eqtl_study_name"_"$gwas_study_name"_coloc_"
	sh run_coloc_for_study_pairing.sh $sig_eqtl_file $all_eqtl_file $gwas_study_file_root $output_root $eqtl_study_name $gwas_study_name "version1"
done<$processed_gwas_studies_file
fi



eqtl_study_name="surge_latent_factor_5_interaction"
sig_eqtl_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_factor_5_interaction_eqtl_results_genome_wide_signficant_bf_fdr_0.05.txt"
all_eqtl_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_factor_5_interaction_eqtl_results_merged.txt"
if false; then
while read gwas_study_name gwas_study_file_root; do
	echo $gwas_study_name
	output_root=$coloc_dir$eqtl_study_name"_"$gwas_study_name"_coloc_"
	sh run_coloc_for_study_pairing.sh $sig_eqtl_file $all_eqtl_file $gwas_study_file_root $output_root $eqtl_study_name $gwas_study_name "version1"
done<$processed_gwas_studies_file
fi


eqtl_study_name="standard_eqtl"
sig_eqtl_file="/work-zfs/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/latent_factor_interaction_eqtl/standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_results_genome_wide_signficant_bf_fdr_0.05.txt"
all_eqtl_file="/work-zfs/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/latent_factor_interaction_eqtl/standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_results_merged.txt"
if false; then
while read gwas_study_name gwas_study_file_root; do
	echo $gwas_study_name
	output_root=$coloc_dir$eqtl_study_name"_"$gwas_study_name"_coloc_"
	sh run_coloc_for_study_pairing.sh $sig_eqtl_file $all_eqtl_file $gwas_study_file_root $output_root $eqtl_study_name $gwas_study_name "version2"
done<$processed_gwas_studies_file
fi


module load R/3.5.1
Rscript visualize_coloc_results.R $processed_gwas_studies_file $coloc_dir $visualize_coloc_dir










