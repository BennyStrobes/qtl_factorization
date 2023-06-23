#!/bin/bash -l
#SBATCH
#SBATCH --time=20:00:00
#SBATCH --mem=20GB
#SBATCH --nodes=1



expression_training_file="$1"
genotype_training_file="$2"
covariate_training_file="$3"
proxy_gene_input_data_dir="$4"
proxy_gene_interaction_eqtl_results_dir="$5"
proxy_gene_jar_file="$6"




#########################
# Convert expression, genotype, and covariates to binary
#########################
# Expression
if false; then
java -jar ${proxy_gene_jar_file} --convertMatrix -i ${expression_training_file} -o ${proxy_gene_input_data_dir}"Expression.binary"
# Genotype
java -jar ${proxy_gene_jar_file} --convertMatrix -i ${genotype_training_file} -o ${proxy_gene_input_data_dir}"Genotypes.binary"
# Covariates
java -jar ${proxy_gene_jar_file} --convertMatrix -i ${covariate_training_file} -o ${proxy_gene_input_data_dir}"Covariates.binary"
fi


#########################
# Run proxy gene analysis
#########################
proxy_gene_input_data_dir2="/scratch16/abattle4/bstrober/qtl_factorization/gtex_cis_eqtl/proxy_gene_input_data"
proxy_gene_interaction_eqtl_results_dir2="/scratch16/abattle4/bstrober/qtl_factorization/gtex_cis_eqtl/proxy_gene_interaction_eqtl"
java -Xmx80g -jar ${proxy_gene_jar_file} -i ${proxy_gene_input_data_dir2} -o ${proxy_gene_interaction_eqtl_results_dir2} -n 1 -pc 5 -ncn -nt "1"
