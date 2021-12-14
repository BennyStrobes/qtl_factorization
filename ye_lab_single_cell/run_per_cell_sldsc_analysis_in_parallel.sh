#!/bin/bash -l

#SBATCH
#SBATCH --time=40:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=15GB





ldsc_source_code_dir="$1"
sldsc_input_data_dir="$2"
sldsc_processed_data_dir="$3"
surge_interaction_eqtl_dir="$4"
sample_names_file="$5"
per_cell_sldsc_processed_data_dir="$6"
per_cell_sldsc_results_dir="$7"
job_number="$8"
num_jobs="$9"

module load python/2.7-anaconda

##########################
# INPUT DATA
###########################
# SURGE loadings across all cells
loading_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_surge_latent_factors_v2.txt"
# SURGE eQTL effect sizes for each latent factor
surge_eqtl_effect_sizes_file=$surge_interaction_eqtl_dir"surge_interaction_eqtl_cis_window_200000_factor_interaction_eqtl_standardized_genotype_results_v2_betas_merged.txt"

##########################
# Run analysis
###########################
# Create sample names file and sample loadings file specific to the parralel run
# This sets up parallelization for the rest of this script
sample_names_file_parr=$per_cell_sldsc_processed_data_dir"sample_names_"$job_number"_"$num_jobs".txt"
loading_file_parr=$per_cell_sldsc_processed_data_dir"surge_sample_loadings_"$job_number"_"$num_jobs".txt"
if false; then
python generate_parrallel_specific_sample_names_and_surge_loadings_for_per_cell_sldsc.py $sample_names_file $loading_file $sample_names_file_parr $loading_file_parr $job_number $num_jobs
fi

# Create sample names file and sample loadings file specific to the parralel run
# This sets up parallelization for the rest of this script
sample_specific_eqtl_effect_sizes_file=$per_cell_sldsc_processed_data_dir"sample_specific_eqtl_effect_sizes_"$job_number"_"$num_jobs".txt"
if false; then
python produce_sample_specific_eqtl_effect_sizes.py $sample_names_file_parr $loading_file_parr $surge_eqtl_effect_sizes_file $sample_specific_eqtl_effect_sizes_file
fi


# Contenate across genes to variant level
# Currently taking sum of squared betas
variant_level_sample_specific_eqtl_effect_sizes_file=$per_cell_sldsc_processed_data_dir"variant_level_sample_specific_eqtl_effect_sizes_"$job_number"_"$num_jobs".txt"
if false; then
python produce_aggregate_sample_specific_eqtl_effect_sizes_to_variant_level.py $sample_specific_eqtl_effect_sizes_file $variant_level_sample_specific_eqtl_effect_sizes_file
fi


# Make annotation file
per_sample_joint_annotation_file_stem=$per_cell_sldsc_processed_data_dir"sample_specific_eqtl_effect_sizes_"$job_number"_"$num_jobs
if false; then
python generate_sample_specific_joint_annotation_file.py $variant_level_sample_specific_eqtl_effect_sizes_file $sldsc_input_data_dir $per_sample_joint_annotation_file_stem
fi

# Generate ld score annotation file
# Ignore ill conditioned error here
if false; then
for chrom_num in $(seq 1 $((22))); do 
	echo $chrom_num
	python ${ldsc_source_code_dir}ldsc.py --l2 --bfile ${sldsc_input_data_dir}"1000G_EUR_Phase3_plink/1000G.EUR.QC."${chrom_num} --ld-wind-cm 1 --annot ${per_sample_joint_annotation_file_stem}"."${chrom_num}".annot" --out ${per_sample_joint_annotation_file_stem}"."${chrom_num} --print-snps ${sldsc_input_data_dir}"hapmap3_snps/hm."${chrom_num}".snp"
done
fi

processed_gwas_studies_file=$sldsc_processed_data_dir"processed_gwas_studies.txt"


if false; then
# Per sample looop
while read sample_name; do
	echo $sample_name
	# Filter annnotation file to eQTL annotations from only this sample
	per_sample_annotation_file_stem=$per_cell_sldsc_processed_data_dir"sample_specific_eqtl_effect_sizes_"$sample_name
	python filter_joint_annotation_files_to_only_one_sample.py $sample_name $sample_names_file_parr $per_sample_joint_annotation_file_stem $per_sample_annotation_file_stem

	# Run sLDSC (on UC) on eQTL annotations from this sample
	gwas_study_name="Ulcerative_Colitis"
	gwas_study_file="/work-zfs/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/sldsc_processed_data/Ulcerative_Colitis.sumstats.gz"
	sldsc_output_stem=$per_cell_sldsc_results_dir"sample_specific_eqtl_effect_sizes_"$gwas_study_name"_"$sample_name
	python ${ldsc_source_code_dir}ldsc.py --h2 ${gwas_study_file} --ref-ld-chr ${sldsc_processed_data_dir}"baselineLD.,"${per_sample_annotation_file_stem}"." --w-ld-chr ${sldsc_input_data_dir}"1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${sldsc_input_data_dir}"1000G_Phase3_frq/1000G.EUR.QC." --out ${sldsc_output_stem}

	# Run sLDSC (on monocyte count) on eQTL annotations from this sample
	gwas_study_name="ukbb_blood_monocyte_count"
	gwas_study_file="/work-zfs/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/sldsc_processed_data/ukbb_blood_monocyte_count.sumstats.gz"
	sldsc_output_stem=$per_cell_sldsc_results_dir"sample_specific_eqtl_effect_sizes_"$gwas_study_name"_"$sample_name
	python ${ldsc_source_code_dir}ldsc.py --h2 ${gwas_study_file} --ref-ld-chr ${sldsc_processed_data_dir}"baselineLD.,"${per_sample_annotation_file_stem}"." --w-ld-chr ${sldsc_input_data_dir}"1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${sldsc_input_data_dir}"1000G_Phase3_frq/1000G.EUR.QC." --out ${sldsc_output_stem}


	# Remove unecessary files
	rm ${per_sample_annotation_file_stem}"."*

done<$sample_names_file_parr
fi














sample_name="1004_1004:0"
per_sample_annotation_file_stem=$per_cell_sldsc_processed_data_dir"sample_specific_eqtl_effect_sizes_"$sample_name
if false; then
python filter_joint_annotation_files_to_only_one_sample.py $sample_name $sample_names_file_parr $per_sample_joint_annotation_file_stem $per_sample_annotation_file_stem
fi
gwas_study_name="Ulcerative_Colitis"
gwas_study_file="/work-zfs/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/sldsc_processed_data/Ulcerative_Colitis.sumstats.gz"
sldsc_output_stem=$per_cell_sldsc_results_dir"sample_specific_eqtl_effect_sizes_"$gwas_study_name"_"$sample_name
python ${ldsc_source_code_dir}ldsc.py --h2 ${gwas_study_file} --ref-ld-chr ${sldsc_processed_data_dir}"baselineLD.,"${per_sample_annotation_file_stem}"." --w-ld-chr ${sldsc_input_data_dir}"1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${sldsc_input_data_dir}"1000G_Phase3_frq/1000G.EUR.QC." --out ${sldsc_output_stem}




