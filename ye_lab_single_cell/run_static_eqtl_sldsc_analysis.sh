#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=lrgmem
#SBATCH --nodes=1
#SBATCH --mem=40GB





ldsc_source_code_dir="$1"
sldsc_input_data_dir="$2"
sldsc_processed_data_dir="$3"
static_eqtl_sldsc_processed_data_dir="$4"
static_eqtl_sldsc_results_dir="$5"
eqtl_effect_sizes_file="$6"


module load python/2.7-anaconda
if false; then
python filter_out_qtl_annotations_from_baselineLD_annotation_file.py ${sldsc_processed_data_dir}"baselineLD." $component_gridspace_sldsc_processed_data_dir"baselineLD_no_qtl."
fi

sample_specific_eqtl_effect_sizes_file=$static_eqtl_sldsc_processed_data_dir"static_eqtl_effect_sizes.txt"
if false; then
python produce_static_eqtl_effect_sizes.py $eqtl_effect_sizes_file $sample_specific_eqtl_effect_sizes_file
fi

# Contenate across genes to variant level
# Currently taking sum of squared betas
variant_level_sample_specific_eqtl_effect_sizes_file=$static_eqtl_sldsc_processed_data_dir"variant_level_static_eqtl_effect_sizes.txt"
if false; then
python produce_aggregate_sample_specific_eqtl_effect_sizes_to_variant_level.py $sample_specific_eqtl_effect_sizes_file $variant_level_sample_specific_eqtl_effect_sizes_file
fi

# Make annotation file
per_sample_joint_annotation_file_stem=$static_eqtl_sldsc_processed_data_dir"static_eqtl_effect_sizes"
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
while read gwas_study_name gwas_study_file; do
	echo $gwas_study_name"_"$gwas_study_file
	sldsc_output_stem=$static_eqtl_sldsc_results_dir"static_eqtl_effect_sizes_"$gwas_study_name
	python ${ldsc_source_code_dir}ldsc.py --h2 ${gwas_study_file} --ref-ld-chr ${static_eqtl_sldsc_processed_data_dir}"baselineLD_no_qtl.,"${static_eqtl_sldsc_processed_data_dir}"static_eqtl_effect_sizes." --w-ld-chr ${sldsc_input_data_dir}"1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." --overlap-annot --print-coefficients --frqfile-chr ${sldsc_input_data_dir}"1000G_Phase3_frq/1000G.EUR.QC." --out ${sldsc_output_stem}
done <$processed_gwas_studies_file
fi





python organize_static_eqtl_sldsc_results_across_parallel_runs.py $static_eqtl_sldsc_processed_data_dir $static_eqtl_sldsc_results_dir







