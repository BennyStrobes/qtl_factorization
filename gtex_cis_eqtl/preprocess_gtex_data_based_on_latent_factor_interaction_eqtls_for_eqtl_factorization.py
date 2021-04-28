import numpy as np 
import os
import sys
import pdb




def extract_eqtl_factorization_tests(cross_tissue_eqtl_results_file, cross_tissue_genome_wide_sig_results_file, num_genes):
	dicti = {}
	binary_arr = []
	temp_dicti = {}
	f = open(cross_tissue_genome_wide_sig_results_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relevent fields
		gene_id = data[1]
		variant_id = data[0]
		test_name = variant_id + ':' + gene_id
		if test_name in dicti:
			print('assumption error')
			pdb.set_trace()
		if len(dicti) < num_genes:
			dicti[test_name] = 1
	f.close()

	f = open(cross_tissue_eqtl_results_file)
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')

		gene_id = data[1]
		variant_id = data[0]
		test_name = variant_id + ':' + gene_id
		pvalue = data[2]
		if test_name in dicti:
			binary_arr.append(1)
			if pvalue != 'NA':
				pvalue = float(pvalue)
				if pvalue > .01:
					print('assumption erroro')
					pdb.set_trace()
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
			t.write(line + '\n')
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
			t.write('\t'.join(standardized_genotype.astype(str)) + '\n')
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


######################
# Command line args
######################
working_dir = sys.argv[1]
num_genes = int(sys.argv[2])


# Input files
latent_factor_interaction_eqtl_genome_wide_sig_results_file = working_dir + 'cross_tissue_latent_factor_interaction_eqtl_results_genome_wide_signficant_bf_fdr_0.05.txt'
latent_factor_interaction_eqtl_results_file = working_dir + 'cross_tissue_latent_factor_interaction_eqtl_results_merged_include_nan.txt'


test_eqtl_dicti, test_eqtl_binary_arr = extract_eqtl_factorization_tests(latent_factor_interaction_eqtl_results_file, latent_factor_interaction_eqtl_genome_wide_sig_results_file, num_genes)


# Generate eqtl factorization test names file
all_test_names_file = working_dir + 'all_tests.txt'
eqtl_factorization_test_names_file = working_dir + 'eqtl_factorization_lf_interaction_eqtl_input_test_names.txt'
generate_eqtl_factorization_test_names_file(all_test_names_file, eqtl_factorization_test_names_file, test_eqtl_binary_arr)

# Generate eqtl factorization gene expression file
#all_gene_expression_file = working_dir + 'cross_tissue_eqtl_residual_expression_input.txt'
all_gene_expression_file = working_dir + 'cross_tissue_eqtl_expression_input.txt'
eqtl_factorization_expression_file = working_dir + 'eqtl_factorization_lf_interaction_eqtl_input_expression.txt'
generate_eqtl_factorization_expression_file(all_gene_expression_file, eqtl_factorization_expression_file, test_eqtl_binary_arr)
save_as_npy_file(eqtl_factorization_expression_file)


# Generate eqtl factorization genotype expression file
all_test_genotype_file = working_dir + 'cross_tissue_eqtl_genotype_input.txt'
eqtl_factorization_genotype_file = working_dir + 'eqtl_factorization_lf_interaction_eqtl_input_genotype.txt'
generate_eqtl_factorization_genotype_file(all_test_genotype_file, eqtl_factorization_genotype_file, test_eqtl_binary_arr)
save_as_npy_file(eqtl_factorization_genotype_file)


# Add intercept to covariate file
covariate_file = working_dir + 'cross_tissue_eqtl_covariate_input.txt'
covariate_with_intercept_file = working_dir + 'cross_tissue_eqtl_covariate_w_intercept_input.txt'
add_intercept_to_covariate_file(covariate_file, covariate_with_intercept_file)
