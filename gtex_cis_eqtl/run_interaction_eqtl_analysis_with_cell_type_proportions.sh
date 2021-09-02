#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1



eqtl_input_dir="$1"
output_stem="$2"


# Input data
qtl_expression_file=$eqtl_input_dir"cross_tissue_eqtl_expression_input.txt"
qtl_genotype_file=$eqtl_input_dir"cross_tissue_eqtl_genotype_input.txt"
qtl_covariate_file=$eqtl_input_dir"cross_tissue_eqtl_covariate_input.txt"
qtl_test_names_file=$eqtl_input_dir"all_tests.txt"
qtl_sample_overlap_file=$eqtl_input_dir"individual_id.txt"
xcell_sample_enrichment_file=$eqtl_input_dir"sample_covariates.txt"


# Loop through cell types
cell_types=( "Adipocytes" "Epithelial_cells" "Hepatocytes" "Keratinocytes" "Myocytes" "Neurons" "Neutrophils")
################################

# Extract xcell environmental variable for each cell type
if false; then
for cell_type in "${cell_types[@]}"; do
	cell_type_context_file=$output_stem"xcell_"$cell_type"_environment_variable.txt"
	python extract_xcell_environmental_variable.py $xcell_sample_enrichment_file $cell_type $cell_type_context_file
done
fi

# Run interaction eqtl analysis in each cell type
num_jobs="1"
job_number="0"
if false; then
for cell_type in "${cell_types[@]}"; do
	cell_type_context_file=$output_stem"xcell_"$cell_type"_environment_variable.txt"
	qtl_output_root=$output_stem$cell_type"_interaction_eqtl_results_"$job_number"_"$num_jobs"_"
	sbatch run_latent_factor_interaction_eqtl_analysis_in_parallel.sh $qtl_expression_file $qtl_genotype_file $qtl_test_names_file $qtl_covariate_file $cell_type_context_file $qtl_sample_overlap_file $qtl_output_root $job_number $num_jobs
done
fi