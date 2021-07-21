import numpy as np 
import os
import sys
import pdb




def extract_eqtl_factorization_tests(test_names_file, num_genes):
	dicti = {}
	binary_arr = []
	temp_dicti = {}

	gene_to_variants = {}
	f = open(test_names_file)
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		variant_id = data[1]
		if gene_id.startswith('HLA'):
			continue
		if gene_id not in gene_to_variants:
			gene_to_variants[gene_id] = []
		gene_to_variants[gene_id].append(variant_id)
	f.close()

	temp_arr = []
	for gene_id in gene_to_variants.keys():
		if len(gene_to_variants[gene_id]) != len(np.unique(gene_to_variants[gene_id])):
			print('assumption erroror')
			continue
		genes_variants = np.asarray(gene_to_variants[gene_id])
		randomly_selected_variant = np.random.choice(genes_variants)
		randomly_selected_test_name = randomly_selected_variant + ';' + gene_id
		temp_arr.append(randomly_selected_test_name)

	randomly_selected_tests = np.random.choice(temp_arr, num_genes, replace=False)
	for test in randomly_selected_tests:
		dicti[test] = 1

	f = open(test_names_file)
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		variant_id = data[1]
		test_name = variant_id + ';' + gene_id
		if test_name in dicti:
			binary_arr.append(1)
		else:
			binary_arr.append(0)
	f.close()
	return dicti, np.asarray(binary_arr)




def generate_eqtl_factorization_test_names_file(all_test_names_file, eqtl_factorization_test_names_file, test_eqtl_binary_arr):
	f = open(all_test_names_file)
	t = open(eqtl_factorization_test_names_file, 'w')
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		if test_eqtl_binary_arr[counter] == 1:
			t.write(line + '\n')
		counter = counter + 1
	f.close()
	t.close()


def generate_eqtl_factorization_expression_file(all_gene_expression_file, eqtl_factorization_expression_file, test_eqtl_binary_arr):
	f = open(all_gene_expression_file)
	t = open(eqtl_factorization_expression_file, 'w')
	counter = 0
	for line in f:
		line = line.rstrip()
		#expr = np.asarray(data).astype(float)
		if test_eqtl_binary_arr[counter] == 1:
			#t.write(line + '\n')
			data = np.asarray(line.split()).astype(float)
			#data[data > 5.0] = 5.0
			#data[data < -5.0] = -5.0
			#data = data - np.mean(data)
			t.write('\t'.join(data.astype(str)) + '\n')
		counter = counter + 1
	t.close()
	f.close()

def generate_eqtl_factorization_genotype_file(all_test_genotype_file, eqtl_factorization_genotype_file, test_eqtl_binary_arr):
	f = open(all_test_genotype_file)
	t = open(eqtl_factorization_genotype_file, 'w')
	counter = 0
	for line in f:
		line = line.rstrip()
		if test_eqtl_binary_arr[counter] == 1:
			data = line.split('\t')
			genotype = np.asarray(data).astype(float)
			standardized_genotype = (genotype - np.mean(genotype))/np.std(genotype)
			t.write('\t'.join(genotype.astype(str)) + '\n')
		counter = counter + 1

	f.close()
	t.close()

def save_as_npy_file(file_name):
	aa = np.loadtxt(file_name)
	npy_file_name = file_name.split('.tx')[0] + '.npy'
	np.save(npy_file_name, aa)

def add_intercept_to_covariate_file(input_file, output_file):
	cov = np.loadtxt(input_file)
	num_samples = cov.shape[0]
	cov_plus_intercept = np.hstack((np.ones((num_samples, 1)), cov))
	np.savetxt(output_file, cov_plus_intercept, fmt="%s", delimiter='\t')


def standardize_columns(cov):
	standardized_cov = np.zeros(cov.shape)
	col_names = np.arange(cov.shape[1])

	for col_name in col_names:
		standardized_cov[:, col_name] = (cov[:, col_name] -  np.mean(cov[:, col_name]))/np.std(cov[:, col_name])
	return standardized_cov

def generate_covariate_file(covariate_file, lf_covariate_file, output_file):
	cov1 = np.loadtxt(covariate_file)
	num_samples = cov1.shape[0]
	cov2 = np.loadtxt(lf_covariate_file)
	num_samples2 = cov2.shape[0]
	if num_samples != num_samples2:
		print('assumption eroorr')
		pdb.set_trace()
	cov = np.hstack((cov1, cov2))
	standardized_cov = standardize_columns(cov)
	cov_plus_intercept = np.hstack((np.ones((num_samples, 1)), standardized_cov))
	np.savetxt(output_file, cov_plus_intercept, fmt="%s", delimiter='\t')


def generate_sample_overlap_file(input_file, output_file):
	f = open(input_file)
	t = open(output_file, 'w')
	for line in f:
		line = line.rstrip()
		t.write(line + '\n')
	f.close()
	t.close()


######################
# Command line args
######################
num_genes = int(sys.argv[1])
working_dir = sys.argv[2]
output_dir = sys.argv[3]



# Input files
#latent_factor_interaction_eqtl_genome_wide_sig_results_file = working_dir + 'standard_eqtl_10.0_qn_zscore_eqtl_results_genome_wide_signficant_bf_fdr_0.2.txt'
#latent_factor_interaction_eqtl_results_file = working_dir + 'latent_factor_interaction_no_cap_5_eqtl_results_merged_include_nan.txt'
latent_factor_interaction_test_names_file = working_dir + 'latent_factor_interaction_10.0_none_zscore_eqtl_input_variant_gene_pairs.txt'

test_eqtl_dicti, test_eqtl_binary_arr = extract_eqtl_factorization_tests(latent_factor_interaction_test_names_file, num_genes)

# Generate eqtl factorization test names file
all_test_names_file = working_dir + 'latent_factor_interaction_10.0_none_zscore_eqtl_input_variant_gene_pairs.txt'
eqtl_factorization_test_names_file = output_dir + 'eqtl_factorization_random_tests_10.0_none_zscore_capped_eqtl_input_test_names.txt'
generate_eqtl_factorization_test_names_file(all_test_names_file, eqtl_factorization_test_names_file, test_eqtl_binary_arr)

# Generate eqtl factorization gene expression file
all_gene_expression_file = working_dir + 'latent_factor_interaction_10.0_none_zscore_eqtl_input_expression.txt'
eqtl_factorization_expression_file = output_dir + 'eqtl_factorization_random_tests_10.0_none_zscore_capped_eqtl_input_expression.txt'
generate_eqtl_factorization_expression_file(all_gene_expression_file, eqtl_factorization_expression_file, test_eqtl_binary_arr)
save_as_npy_file(eqtl_factorization_expression_file)


# Generate eqtl factorization genotype expression file
all_test_genotype_file = working_dir + 'latent_factor_interaction_10.0_none_zscore_eqtl_input_genotype.txt'
eqtl_factorization_genotype_file = output_dir + 'eqtl_factorization_random_tests_10.0_none_zscore_capped_eqtl_input_unnormalized_genotype.txt'
generate_eqtl_factorization_genotype_file(all_test_genotype_file, eqtl_factorization_genotype_file, test_eqtl_binary_arr)
save_as_npy_file(eqtl_factorization_genotype_file)

# Generate covariate file
covariate_file = working_dir + 'latent_factor_interaction_10.0_none_zscore_eqtl_input_covariates.txt'
lf_covariate_file = working_dir + 'latent_factor_interaction_10.0_none_zscore_eqtl_input_latent_factors.txt'
covariate_with_intercept_file = output_dir + 'eqtl_factorization_random_tests_10.0_none_zscore_capped_eqtl_input_covariates.txt'
generate_covariate_file(covariate_file, lf_covariate_file, covariate_with_intercept_file)

# Generate sample overlap file
sample_overlap_input_file = working_dir + 'latent_factor_interaction_10.0_none_zscore_eqtl_input_sample_overlap.txt'
sample_overlap_output_file = output_dir + 'eqtl_factorization_random_tests_10.0_none_zscore_capped_eqtl_input_sample_overlap.txt'
generate_sample_overlap_file(sample_overlap_input_file, sample_overlap_output_file)




'''
# Input files
latent_factor_interaction_eqtl_genome_wide_sig_results_file = working_dir + 'latent_factor_interaction_10.0_qn_zscore_eqtl_results_genome_wide_signficant_bf_fdr_0.05.txt'
#latent_factor_interaction_eqtl_results_file = working_dir + 'latent_factor_interaction_no_cap_5_eqtl_results_merged_include_nan.txt'
latent_factor_interaction_test_names_file = working_dir + 'latent_factor_interaction_10.0_qn_zscore_eqtl_input_variant_gene_pairs.txt'

test_eqtl_dicti, test_eqtl_binary_arr = extract_eqtl_factorization_tests(latent_factor_interaction_test_names_file, latent_factor_interaction_eqtl_genome_wide_sig_results_file, num_genes)

# Generate eqtl factorization test names file
all_test_names_file = working_dir + 'latent_factor_interaction_10.0_qn_zscore_eqtl_input_variant_gene_pairs.txt'
eqtl_factorization_test_names_file = output_dir + 'eqtl_factorization_lf_interaction_10.0_qn_zscore_capped_eqtl_input_test_names.txt'
generate_eqtl_factorization_test_names_file(all_test_names_file, eqtl_factorization_test_names_file, test_eqtl_binary_arr)

# Generate eqtl factorization gene expression file
all_gene_expression_file = working_dir + 'latent_factor_interaction_10.0_qn_zscore_eqtl_input_expression.txt'
eqtl_factorization_expression_file = output_dir + 'eqtl_factorization_lf_interaction_10.0_qn_zscore_capped_eqtl_input_expression.txt'
generate_eqtl_factorization_expression_file(all_gene_expression_file, eqtl_factorization_expression_file, test_eqtl_binary_arr)
save_as_npy_file(eqtl_factorization_expression_file)


# Generate eqtl factorization genotype expression file
all_test_genotype_file = working_dir + 'latent_factor_interaction_10.0_qn_zscore_eqtl_input_genotype.txt'
eqtl_factorization_genotype_file = output_dir + 'eqtl_factorization_lf_interaction_10.0_qn_zscore_capped_eqtl_input_genotype.txt'
generate_eqtl_factorization_genotype_file(all_test_genotype_file, eqtl_factorization_genotype_file, test_eqtl_binary_arr)
save_as_npy_file(eqtl_factorization_genotype_file)

# Generate covariate file
covariate_file = working_dir + 'latent_factor_interaction_10.0_qn_zscore_eqtl_input_covariates.txt'
lf_covariate_file = working_dir + 'latent_factor_interaction_10.0_qn_zscore_eqtl_input_latent_factors.txt'
covariate_with_intercept_file = output_dir + 'eqtl_factorization_lf_interaction_10.0_qn_zscore_capped_eqtl_input_covariates.txt'
generate_covariate_file(covariate_file, lf_covariate_file, covariate_with_intercept_file)

# Generate sample overlap file
sample_overlap_input_file = working_dir + 'latent_factor_interaction_10.0_qn_zscore_eqtl_input_sample_overlap.txt'
sample_overlap_output_file = output_dir + 'eqtl_factorization_lf_interaction_10.0_qn_zscore_capped_eqtl_input_sample_overlap.txt'
generate_sample_overlap_file(sample_overlap_input_file, sample_overlap_output_file)

'''

