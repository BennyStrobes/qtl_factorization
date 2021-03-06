#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --mem=17GB





ldsc_source_code_dir="$1"
sldsc_input_data_dir="$2"
sldsc_processed_data_dir="$3"
surge_interaction_eqtl_dir="$4"
sample_names_file="$5"
per_cell_sldsc_processed_data_dir="$6"
per_cell_sldsc_results_dir="$7"
custom_ldsc_source_code_dir="$8"
loading_file="$9"
surge_eqtl_effect_sizes_file="${10}"


module load python/2.7-anaconda

if false; then
python filter_out_qtl_annotations_from_baselineLD_annotation_file.py ${sldsc_processed_data_dir}"baselineLD." $per_cell_sldsc_processed_data_dir"baselineLD_no_qtl."
fi




if false; then
job_number="0"
num_jobs="100"
sbatch run_per_cell_sldsc_analysis_in_parallel.sh $ldsc_source_code_dir $sldsc_input_data_dir $sldsc_processed_data_dir $surge_interaction_eqtl_dir $sample_names_file $per_cell_sldsc_processed_data_dir $per_cell_sldsc_results_dir $job_number $num_jobs $custom_ldsc_source_code_dir $loading_file $surge_eqtl_effect_sizes_file

job_number="1"
num_jobs="100"
sbatch run_per_cell_sldsc_analysis_in_parallel.sh $ldsc_source_code_dir $sldsc_input_data_dir $sldsc_processed_data_dir $surge_interaction_eqtl_dir $sample_names_file $per_cell_sldsc_processed_data_dir $per_cell_sldsc_results_dir $job_number $num_jobs $custom_ldsc_source_code_dir $loading_file $surge_eqtl_effect_sizes_file

job_number="2"
num_jobs="100"
sbatch run_per_cell_sldsc_analysis_in_parallel.sh $ldsc_source_code_dir $sldsc_input_data_dir $sldsc_processed_data_dir $surge_interaction_eqtl_dir $sample_names_file $per_cell_sldsc_processed_data_dir $per_cell_sldsc_results_dir $job_number $num_jobs $custom_ldsc_source_code_dir $loading_file $surge_eqtl_effect_sizes_file

job_number="3"
num_jobs="100"
sbatch run_per_cell_sldsc_analysis_in_parallel.sh $ldsc_source_code_dir $sldsc_input_data_dir $sldsc_processed_data_dir $surge_interaction_eqtl_dir $sample_names_file $per_cell_sldsc_processed_data_dir $per_cell_sldsc_results_dir $job_number $num_jobs $custom_ldsc_source_code_dir $loading_file $surge_eqtl_effect_sizes_file

job_number="4"
num_jobs="100"
sbatch run_per_cell_sldsc_analysis_in_parallel.sh $ldsc_source_code_dir $sldsc_input_data_dir $sldsc_processed_data_dir $surge_interaction_eqtl_dir $sample_names_file $per_cell_sldsc_processed_data_dir $per_cell_sldsc_results_dir $job_number $num_jobs $custom_ldsc_source_code_dir $loading_file $surge_eqtl_effect_sizes_file

job_number="5"
num_jobs="100"
sbatch run_per_cell_sldsc_analysis_in_parallel.sh $ldsc_source_code_dir $sldsc_input_data_dir $sldsc_processed_data_dir $surge_interaction_eqtl_dir $sample_names_file $per_cell_sldsc_processed_data_dir $per_cell_sldsc_results_dir $job_number $num_jobs $custom_ldsc_source_code_dir $loading_file $surge_eqtl_effect_sizes_file
fi



num_jobs="100"
if false; then
for job_number in $(seq 0 $(($num_jobs-1))); do 
	sbatch run_per_cell_sldsc_analysis_in_parallel.sh $ldsc_source_code_dir $sldsc_input_data_dir $sldsc_processed_data_dir $surge_interaction_eqtl_dir $sample_names_file $per_cell_sldsc_processed_data_dir $per_cell_sldsc_results_dir $job_number $num_jobs $custom_ldsc_source_code_dir $loading_file $surge_eqtl_effect_sizes_file
done
fi


if false; then
python organize_per_cell_sldsc_results_across_parallel_runs.py $sample_names_file $per_cell_sldsc_results_dir
fi

