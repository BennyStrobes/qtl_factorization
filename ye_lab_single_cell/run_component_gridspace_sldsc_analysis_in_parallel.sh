#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --mem=40GB





ldsc_source_code_dir="$1"
sldsc_input_data_dir="$2"
sldsc_processed_data_dir="$3"
sample_names_file="$4"
per_cell_sldsc_processed_data_dir="$5"
per_cell_sldsc_results_dir="$6"
job_number="$7"
num_jobs="$8"
custom_ldsc_source_code_dir="${9}"
loading_file="${10}"
surge_eqtl_effect_sizes_file="${11}"


source ~/.bash_profile
conda activate ldsc



##########################
# Run analysis
###########################
# Create sample names file specific to the parralel run
# This sets up parallelization for the rest of this script
sample_names_file_parr=$per_cell_sldsc_processed_data_dir"sample_names_"$job_number"_"$num_jobs".txt"
python generate_parrallel_specific_sample_names_for_component_gridspace_sldsc.py $sample_names_file $sample_names_file_parr $job_number $num_jobs


# Create sample names file and sample loadings file specific to the parralel run
# This sets up parallelization for the rest of this script
sample_specific_eqtl_effect_sizes_file=$per_cell_sldsc_processed_data_dir"sample_specific_eqtl_effect_sizes_"$job_number"_"$num_jobs".txt"
python produce_component_gridspace_sample_specific_eqtl_effect_sizes.py $sample_names_file_parr $surge_eqtl_effect_sizes_file $sample_specific_eqtl_effect_sizes_file


# Contenate across genes to variant level
# Currently taking sum of squared betas
variant_level_sample_specific_eqtl_effect_sizes_file=$per_cell_sldsc_processed_data_dir"variant_level_sample_specific_eqtl_effect_sizes_"$job_number"_"$num_jobs".txt"
python produce_aggregate_sample_specific_eqtl_effect_sizes_to_variant_level.py $sample_specific_eqtl_effect_sizes_file $variant_level_sample_specific_eqtl_effect_sizes_file


# Make annotation file
per_sample_joint_annotation_file_stem=$per_cell_sldsc_processed_data_dir"sample_specific_eqtl_effect_sizes_"$job_number"_"$num_jobs
python generate_sample_specific_joint_annotation_file.py $variant_level_sample_specific_eqtl_effect_sizes_file $sldsc_input_data_dir $per_sample_joint_annotation_file_stem


# Generate ld score annotation file
# Ignore ill conditioned error here (because not actually going to model all annotations jointly)
for chrom_num in $(seq 1 $((22))); do 
	echo $chrom_num
	python ${ldsc_source_code_dir}ldsc.py --l2 --bfile ${sldsc_input_data_dir}"1000G_EUR_Phase3_plink/1000G.EUR.QC."${chrom_num} --ld-wind-cm 1 --annot ${per_sample_joint_annotation_file_stem}"."${chrom_num}".annot" --out ${per_sample_joint_annotation_file_stem}"."${chrom_num} --print-snps ${sldsc_input_data_dir}"hapmap3_snps/hm."${chrom_num}".snp"
done



# Generate sample specific .npy annotation files
per_sample_annotation_summary_file=$per_cell_sldsc_processed_data_dir"summary_sample_specific_eqtl_effect_sizes_"$job_number"_"$num_jobs".txt"
python filter_joint_annotation_files_to_only_one_file_per_sample.py $sample_names_file_parr $per_sample_joint_annotation_file_stem $per_sample_annotation_summary_file



#processed_gwas_studies_file=$sldsc_processed_data_dir"processed_gwas_studies_short.txt"
processed_gwas_studies_file=$sldsc_processed_data_dir"processed_gwas_studies.txt"
while read gwas_study_name gwas_study_file; do
	echo $gwas_study_name"_"$gwas_study_file
	sldsc_output_stem=$per_cell_sldsc_results_dir"sample_specific_eqtl_effect_sizes_"$gwas_study_name
	sldsc_log_output_stem=$per_cell_sldsc_results_dir"sample_specific_eqtl_effect_sizes_"$gwas_study_name"_"$job_number"_"$num_jobs

	python ${custom_ldsc_source_code_dir}ldsc.py --h2 ${gwas_study_file} --ref-ld-chr ${per_cell_sldsc_processed_data_dir}"baselineLD_no_qtl." --multi_single_anno $per_sample_annotation_summary_file --w-ld-chr ${sldsc_input_data_dir}"1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${sldsc_input_data_dir}"1000G_Phase3_frq/1000G.EUR.QC." --out ${sldsc_output_stem} --logout ${sldsc_log_output_stem}
done <$processed_gwas_studies_file




