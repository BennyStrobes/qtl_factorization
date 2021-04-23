


#########################
# Output directories
#########################
# Root directory for simulation analysis
root_dir="/work-zfs/abattle4/bstrober/qtl_factorization/gtex_cis_eqtl/"
# Directory containing input data
input_data_dir=$root_dir"input_data/"
# Directory containing simulated data
processed_data_dir=$root_dir"processed_data/"
# Directory containing eqtl results on simulated data
eqtl_results_dir=$root_dir"eqtl_results/"
# Directory containing visualizations
visualization_expression_dir=$root_dir"visualize_expression/"
# Directory containing visualizations of standard eqtl approachees
visualize_standard_eqtl_dir=$root_dir"visualize_standard_eqtl/"
# Directory containing visualization of eQTL Factorization results
visualize_eqtl_factorization_results_dir=$root_dir"visualize_eqtl_factorization/"

#########################
# Input Data
#########################
gtex_expression_dir="/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_expression_matrices/"
gtex_tpm_dir="/work-zfs/abattle4/lab_data/GTEx_v8/processed/rna_seq_by_tissue/gene_tpm/"
gtex_covariate_dir="/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_covariates/"
gtex_genotype_dir="/work-zfs/abattle4/lab_data/GTEx_v8_trans_eqtl_data_processed_by_brian/processed_genotypes/"
gtex_egene_dir="/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL/"
gtex_eqtl_dir="/work-zfs/abattle4/lab_data/GTEx_v8/ciseQTL/GTEx_Analysis_v8_eQTL_all_associations/"
gtex_tissue_colors_file="/work-zfs/abattle4/bstrober/single_cell_eqtl_factorization/gtex_cis_eqtl/input_data/gtex_colors.txt"
gtex_individual_information_file="/work-zfs/abattle4/lab_data/GTEx_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt"
gtex_sample_information_file="/work-zfs/abattle4/bstrober/qtl_factorization/gtex_cis_eqtl/input_data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
cell_type_decomposition_hlv_file="/work-zfs/abattle4/marios/GTEx_v8/coloc/scRNA_decon/CIBERSORT_GTEx_LV.csv"
african_ancestry_gtex_eqtl_dir="/work-zfs/abattle4/surya/worksets/for_ben/GTEx_Analysis_v8_AFR_eQTL/"
european_ancestry_gtex_eqtl_dir="/work-zfs/abattle4/surya/worksets/for_ben/GTEx_Analysis_v8_EUR_eQTL/"
gtex_xcell_enrichment_file="/work-zfs/abattle4/lab_data/GTEx_v8_cell_type_interaction/GTEx_Analysis_v8_xCell_scores_7_celltypes.txt"







######################
## 1 tissue (whole blood) case
########################
tissues_file=$input_data_dir"tissues_subset_whole_blood.txt"
output_dir=$processed_data_dir"tissues_subset_whole_blood_"
output_visualization_dir=$visualization_expression_dir"tissues_subset_whole_blood_"
output_eqtl_visualization_dir=$visualize_standard_eqtl_dir"tissues_subset_whole_blood_"
sh preprocess_gtex_data_for_eqtl_factorization.sh $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $gtex_individual_information_file $gtex_sample_information_file $gtex_eqtl_dir $gtex_xcell_enrichment_file $output_dir $output_visualization_dir $output_eqtl_visualization_dir




## Run eqtl factorization!
sample_overlap_file=$processed_data_dir"tissues_subset_10_individual_id.txt"
expression_training_file=$processed_data_dir"tissues_subset_10_eqtl_factorization_lf_interaction_eqtl_input_expression.npy"
genotype_training_file=$processed_data_dir"tissues_subset_10_eqtl_factorization_lf_interaction_eqtl_input_genotype.npy"
covariate_file=$processed_data_dir"tissues_subset_10_cross_tissue_eqtl_residual_covariate_w_intercept_input.txt"
num_latent_factors="20"
lambda_v="1"
model_name="eqtl_factorization_vi"
seed="6"

output_stem=$eqtl_results_dir"tissues_subset_10_lf_interaction_egenes_"$model_name"_results_k_init_"$num_latent_factors"_lambda_v_"$lambda_v"_seed_"$seed"_RE_init2_"
if false; then
sh run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem
fi

if false; then
module load R/3.5.1
Rscript visualize_eqtl_factorization.R $processed_data_dir $eqtl_results_dir $visualize_eqtl_factorization_results_dir $gtex_tissue_colors_file
fi









#####################
## 10 tissues case
#####################
tissues_file=$input_data_dir"tissues_subset_10.txt"
output_dir=$processed_data_dir"tissues_subset_10_"
output_visualization_dir=$visualization_expression_dir"tissues_subset_10_"
output_eqtl_visualization_dir=$visualize_standard_eqtl_dir"tissues_subset_10_"
if false; then
sh preprocess_gtex_data_for_eqtl_factorization.sh $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $gtex_individual_information_file $gtex_sample_information_file $gtex_eqtl_dir $gtex_xcell_enrichment_file $output_dir $output_visualization_dir $output_eqtl_visualization_dir
fi



## Run eqtl factorization!
sample_overlap_file=$processed_data_dir"tissues_subset_10_individual_id.txt"
expression_training_file=$processed_data_dir"tissues_subset_10_eqtl_factorization_lf_interaction_eqtl_input_expression.npy"
genotype_training_file=$processed_data_dir"tissues_subset_10_eqtl_factorization_lf_interaction_eqtl_input_genotype.npy"
covariate_file=$processed_data_dir"tissues_subset_10_cross_tissue_eqtl_residual_covariate_w_intercept_input.txt"
num_latent_factors="20"
lambda_v="1"
model_name="eqtl_factorization_vi"
seed="6"

output_stem=$eqtl_results_dir"tissues_subset_10_lf_interaction_egenes_"$model_name"_results_k_init_"$num_latent_factors"_lambda_v_"$lambda_v"_seed_"$seed"_RE_init2_"
if false; then
sh run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $seed $output_stem
fi

if false; then
module load R/3.5.1
Rscript visualize_eqtl_factorization.R $processed_data_dir $eqtl_results_dir $visualize_eqtl_factorization_results_dir $gtex_tissue_colors_file
fi
