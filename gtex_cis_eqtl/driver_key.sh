


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
visualization_dir=$root_dir"visualize/"



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
cell_type_decomposition_hlv_file="/work-zfs/abattle4/marios/GTEx_v8/coloc/scRNA_decon/CIBERSORT_GTEx_LV.csv"
african_ancestry_gtex_eqtl_dir="/work-zfs/abattle4/surya/worksets/for_ben/GTEx_Analysis_v8_AFR_eQTL/"
european_ancestry_gtex_eqtl_dir="/work-zfs/abattle4/surya/worksets/for_ben/GTEx_Analysis_v8_EUR_eQTL/"



## 10 tissues case
tissues_file=$input_data_dir"tissues_subset_10.txt"
output_dir=$processed_data_dir"tissues_subset_10_"
if false; then
sh preprocess_gtex_data_for_eqtl_factorization.sh $tissues_file $gtex_expression_dir $gtex_tpm_dir $gtex_covariate_dir $gtex_genotype_dir $gtex_egene_dir $gtex_individual_information_file $gtex_eqtl_dir $output_dir
fi




## Run eqtl factorization!
sample_overlap_file=$processed_data_dir"tissues_subset_10_individual_id.txt"
expression_training_file=$processed_data_dir"tissues_subset_10_eqtl_factorization_p_thresh_0.01_input_expression.npy"
genotype_training_file=$processed_data_dir"tissues_subset_10_eqtl_factorization_p_thresh_0.01_input_genotype.npy"
covariate_file=$processed_data_dir"tissues_subset_10_cross_tissue_eqtl_covariate_w_intercept_input.txt"
num_latent_factors="10"
lambda_v="1"
model_name="eqtl_factorization_vi"


output_stem=$eqtl_results_dir"tissues_subset_10_gene_thresh_0.01_"$model_name"_results_k_init_"$num_latent_factors"_lambda_v_"$lambda_v"_"
sh run_eqtl_factorization.sh $expression_training_file $genotype_training_file $covariate_file $sample_overlap_file $num_latent_factors $lambda_v $model_name $output_stem

