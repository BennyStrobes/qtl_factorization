#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=20GB

sldsc_input_data_dir="$1"
surge_interaction_eqtl_dir="$2"
sldsc_processed_data_dir="$3"
ldsc_source_code_dir="$4"

if false; then
python generate_sldsc_annotation_file_based_on_all_surge_eqtls.py $sldsc_input_data_dir $surge_interaction_eqtl_dir $sldsc_processed_data_dir
fi

if false; then
for component_num in $(seq 1 $((10))); do 
	echo $component_num
	python generate_sldsc_annotation_file_based_on_all_single_component_surge_eqtls.py $sldsc_input_data_dir $surge_interaction_eqtl_dir $sldsc_processed_data_dir $component_num
done
fi

if false; then
python generate_sldsc_annotation_file_based_on_all_standard_eqtls.py $sldsc_input_data_dir $surge_interaction_eqtl_dir $sldsc_processed_data_dir
fi



if false; then
python generate_sldsc_annotation_file_based_on_surge_eqtls.py $sldsc_input_data_dir $surge_interaction_eqtl_dir $sldsc_processed_data_dir

python generate_sldsc_annotation_file_based_on_surge_eqtls_and_base.py $sldsc_input_data_dir $surge_interaction_eqtl_dir $sldsc_processed_data_dir
fi
if false; then
python generate_sldsc_annotation_file_based_on_real_valued_surge_eqtls.py $sldsc_input_data_dir $surge_interaction_eqtl_dir $sldsc_processed_data_dir
fi


if false; then
for chrom_num in $(seq 1 $((22))); do 
	sbatch generate_annotation_weighted_ld_scores_in_chromosome.sh $ldsc_source_code_dir $sldsc_input_data_dir $sldsc_processed_data_dir $chrom_num
done 
fi






