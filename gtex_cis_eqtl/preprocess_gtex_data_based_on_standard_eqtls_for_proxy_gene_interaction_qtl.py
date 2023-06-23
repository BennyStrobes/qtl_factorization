import numpy as np 
import os
import sys
import pdb




def extract_eqtl_factorization_tests(cross_tissue_eqtl_results_file, cross_tissue_genome_wide_sig_results_file):
	dicti = {}
	binary_arr = []
	temp_dicti = {}
	f = open(cross_tissue_genome_wide_sig_results_file)
	head_count = 0
	used_genes = {}
	used_variants = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Extract relevent fields
		gene_id = data[1]
		variant_id = data[0]
		if gene_id in used_genes:
			continue
		if variant_id in used_variants:
			continue
		test_name = variant_id + ':' + gene_id
		if test_name in dicti:
			print('assumption error')
			pdb.set_trace()
		dicti[test_name] = 1
		used_genes[gene_id] = 1
		used_variants[variant_id] = 1
	f.close()

	f = open(cross_tissue_eqtl_results_file)
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')

		gene_id = data[0]
		variant_id = data[1]
		test_name = variant_id + ':' + gene_id
		if head_count == 0:
			head_count = head_count + 1
			continue
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


def generate_eqtl_factorization_expression_file(all_gene_expression_file, eqtl_factorization_expression_file, test_eqtl_binary_arr, gene_names, sample_names):
	f = open(all_gene_expression_file)
	t = open(eqtl_factorization_expression_file, 'w')
	t.write('-\t' + '\t'.join(sample_names) + '\n')
	counter = 0
	gene_counter = 0
	for line in f:
		line = line.rstrip()
		#data = line.split('\t')
		#expr = np.asarray(data).astype(float)
		if test_eqtl_binary_arr[counter] == 1:
			t.write(gene_names[gene_counter] + '\t' + line + '\n')
			gene_counter = gene_counter + 1
		counter = counter + 1
	t.close()
	f.close()
	if len(gene_names) != gene_counter:
		print('assumption erororro')
		pdb.set_trace()
	return

def generate_eqtl_factorization_genotype_file(all_test_genotype_file, eqtl_factorization_genotype_file, test_eqtl_binary_arr, variant_names, sample_names):
	f = open(all_test_genotype_file)
	t = open(eqtl_factorization_genotype_file, 'w')
	t.write('-\t' + '\t'.join(sample_names) + '\n')
	counter = 0
	variant_counter = 0
	for line in f:
		line = line.rstrip()
		if test_eqtl_binary_arr[counter] == 1:
			#data = line.split('\t')
			#genotype = np.asarray(data)
			#standardized_genotype = (genotype - np.mean(genotype))/np.std(genotype)
			t.write(variant_names[variant_counter] + '\t' + line + '\n')
			variant_counter = variant_counter + 1
		counter = counter + 1

	if len(variant_names) != variant_counter:
		print('assumption erororo')
		pdb.set_trace()

	f.close()
	t.close()
	return

def save_as_npy_file(file_name):
	aa = np.loadtxt(file_name)
	npy_file_name = file_name.split('.tx')[0] + '.npy'
	np.save(npy_file_name, aa)

def add_intercept_to_covariate_file(input_file, output_file):
	cov = np.loadtxt(input_file)
	num_samples = cov.shape[0]
	cov_plus_intercept = np.hstack((np.ones((num_samples, 1)), cov))
	np.savetxt(output_file, cov_plus_intercept, fmt="%s", delimiter='\t')

def extract_sample_names(filer):
	f = open(filer)
	arr = []
	for line in f:
		line = line.rstrip()
		arr.append(line)
	f.close()
	return np.asarray(arr)

def extract_gene_names_and_variant_names_from_test_names_file(filer):
	f = open(filer)
	genes = []
	variants = []
	tests = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		if head_count == 0:
			head_count = head_count + 1
			continue
		data = line.split('\t')
		gene_name = data[0]
		variant_name = data[1]
		genes.append(gene_name)
		variants.append(variant_name)
	f.close()

	return np.asarray(genes), np.asarray(variants)

def get_covariate_names(filer):
	f = open(filer)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			covs = data[1:]
			continue
	f.close()

	return covs

def create_proxy_gene_covariate_file(covariate_file, new_covariate_file, covariate_names, sample_names, raw_expression_file):
	old_cov = np.transpose(np.loadtxt(covariate_file,dtype=str,delimiter='\t'))
	t = open(new_covariate_file,'w')
	t.write('-\t' + '\t'.join(sample_names) + '\n')

	if len(covariate_names) != old_cov.shape[0]:
		print('assumption eroror')
		pdb.set_trace()

	for cov_iter, cov_name in enumerate(covariate_names):
		t.write(cov_name + '\t' + '\t'.join(old_cov[cov_iter, :]) + '\n')
	
	head_count = 0
	f = open(raw_expression_file)
	for line in f:
		line = line.rstrip()
		if head_count == 0:
			data = line.split('\t')
			head_count = head_count + 1
			line_sample_names = data[1:]
			if np.array_equal(line_sample_names, sample_names) == False:
				print('asssumption eroror')
				pdb.set_trace()
			continue
		t.write(line + '\n')
	f.close()
	t.close()
	return


######################
# Command line args
######################
working_dir = sys.argv[1]
qtl_model_version = sys.argv[2]

# Sample names
sample_names = extract_sample_names(working_dir+'outliers_removed_sample_names.txt')

# Input files
latent_factor_interaction_eqtl_genome_wide_sig_results_file = working_dir + 'cross_tissue_eqtl_results_genome_wide_signficant_bf_fdr_0.05_' + qtl_model_version + '_results.txt'
latent_factor_interaction_eqtl_results_file = working_dir + 'all_tests.txt'


test_eqtl_dicti, test_eqtl_binary_arr = extract_eqtl_factorization_tests(latent_factor_interaction_eqtl_results_file, latent_factor_interaction_eqtl_genome_wide_sig_results_file)

# Generate eqtl factorization test names file
all_test_names_file = working_dir + 'all_tests.txt'
eqtl_factorization_test_names_file = working_dir + 'proxy_gene_factorization_standard_eqtl_' + qtl_model_version + '_input_test_names.txt'
generate_eqtl_factorization_test_names_file(all_test_names_file, eqtl_factorization_test_names_file, test_eqtl_binary_arr)

# Get gene names and variant names
gene_names, variant_names = extract_gene_names_and_variant_names_from_test_names_file(eqtl_factorization_test_names_file)


# Generate eqtl factorization gene expression file
#all_gene_expression_file = working_dir + 'cross_tissue_eqtl_residual_expression_input.txt'
all_gene_expression_file = working_dir + 'cross_tissue_eqtl_expression_input.txt'
eqtl_factorization_expression_file = working_dir + 'proxy_gene_factorization_standard_eqtl_' + qtl_model_version + '_input_expression.txt'
generate_eqtl_factorization_expression_file(all_gene_expression_file, eqtl_factorization_expression_file, test_eqtl_binary_arr, gene_names, sample_names)

#print('genotype')
# Generate eqtl factorization genotype expression file
all_test_genotype_file = working_dir + 'cross_tissue_eqtl_genotype_input.txt'
eqtl_factorization_genotype_file = working_dir + 'proxy_gene_factorization_standard_eqtl_' + qtl_model_version + '_input_genotype.txt'
generate_eqtl_factorization_genotype_file(all_test_genotype_file, eqtl_factorization_genotype_file, test_eqtl_binary_arr, variant_names, sample_names)

# Raw expression
raw_expression_file = working_dir + 'outliers_removed_cross_tissue_tpm_standardized.txt'


#print('covariate')
# Add intercept to covariate file
covariate_names = get_covariate_names(working_dir + 'covariates.txt')
covariate_file = working_dir + 'cross_tissue_eqtl_covariate_input.txt'
new_covariate_file = working_dir + 'proxy_gene_factorization_' + qtl_model_version + '_covariate_input.txt'
create_proxy_gene_covariate_file(covariate_file, new_covariate_file, covariate_names, sample_names, raw_expression_file)



