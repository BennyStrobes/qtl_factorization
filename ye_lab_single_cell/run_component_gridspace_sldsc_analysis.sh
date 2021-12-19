#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=17GB





ldsc_source_code_dir="$1"
sldsc_input_data_dir="$2"
sldsc_processed_data_dir="$3"
component_gridspace_sldsc_processed_data_dir="$4"
component_gridspace_sldsc_results_dir="$5"
custom_ldsc_source_code_dir="$6"
loading_file="$7"
surge_eqtl_effect_sizes_file="$8"


module load python/2.7-anaconda
if false; then
python filter_out_qtl_annotations_from_baselineLD_annotation_file.py ${sldsc_processed_data_dir}"baselineLD." $component_gridspace_sldsc_processed_data_dir"baselineLD_no_qtl."
fi


num_samples="200"
sample_names_file=$component_gridspace_sldsc_processed_data_dir"component_gridspace_"$num_samples"_sample_names.txt"
if false; then
python create_component_gridspace_sample_names.py $loading_file $sample_names_file $num_samples
fi


num_jobs="5"
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	sbatch run_component_gridspace_sldsc_analysis_in_parallel.sh $ldsc_source_code_dir $sldsc_input_data_dir $sldsc_processed_data_dir $sample_names_file $component_gridspace_sldsc_processed_data_dir $component_gridspace_sldsc_results_dir $job_number $num_jobs $custom_ldsc_source_code_dir $loading_file $surge_eqtl_effect_sizes_file
done
fi





python organize_component_gridspace_sldsc_results_across_parallel_runs.py $sample_names_file $component_gridspace_sldsc_results_dir







