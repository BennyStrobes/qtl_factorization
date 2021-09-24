#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1




surge_interaction_eqtl_dir="$1"
gwas_overlap_dir="$2"
gwas_files="$3"
test_info_file="$4"
gwas_snp_info_dir="$5"




# Standard eqtls
if false; then

egene_file="/work-zfs/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/latent_factor_interaction_eqtl/standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_eqtl_results_genome_wide_signficant_bf_fdr_0.05.txt"
output_root=$gwas_overlap_dir"gwas_overlap_standard_eqtls_"
echo $factor_num
python check_gwas_overlap.py $egene_file $gwas_files $gwas_snp_info_dir $output_root
fi

factor_num="1"
egene_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_factor_"$factor_num"_interaction_eqtl_results_genome_wide_signficant_bf_fdr_0.1.txt"
output_root=$gwas_overlap_dir"gwas_overlap_surge_context_"$factor_num"_"
python check_gwas_overlap.py $egene_file $gwas_files $test_info_file $gwas_snp_info_dir $output_root



if false; then

factor_num="2"
egene_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_factor_"$factor_num"_interaction_eqtl_results_genome_wide_signficant_bf_fdr_0.1.txt"
output_root=$gwas_overlap_dir"gwas_overlap_surge_context_"$factor_num"_"
python check_gwas_overlap3.py $egene_file $gwas_files $test_info_file $gwas_snp_info_dir $output_root

factor_num="3"
egene_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_factor_"$factor_num"_interaction_eqtl_results_genome_wide_signficant_bf_fdr_0.1.txt"
output_root=$gwas_overlap_dir"gwas_overlap_surge_context_"$factor_num"_"
python check_gwas_overlap3.py $egene_file $gwas_files $test_info_file $gwas_snp_info_dir $output_root

factor_num="4"
egene_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_factor_"$factor_num"_interaction_eqtl_results_genome_wide_signficant_bf_fdr_0.1.txt"
output_root=$gwas_overlap_dir"gwas_overlap_surge_context_"$factor_num"_"
python check_gwas_overlap3.py $egene_file $gwas_files $test_info_file $gwas_snp_info_dir $output_root
factor_num="6"
egene_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_factor_"$factor_num"_interaction_eqtl_results_genome_wide_signficant_bf_fdr_0.1.txt"
output_root=$gwas_overlap_dir"gwas_overlap_surge_context_"$factor_num"_"
python check_gwas_overlap3.py $egene_file $gwas_files $test_info_file $gwas_snp_info_dir $output_root

factor_num="7"
egene_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_factor_"$factor_num"_interaction_eqtl_results_genome_wide_signficant_bf_fdr_0.1.txt"
output_root=$gwas_overlap_dir"gwas_overlap_surge_context_"$factor_num"_"
python check_gwas_overlap3.py $egene_file $gwas_files $test_info_file $gwas_snp_info_dir $output_root

fi




if false; then
module load R/3.5.1
Rscript visualize_gwas_overlap.R $gwas_overlap_dir
fi

