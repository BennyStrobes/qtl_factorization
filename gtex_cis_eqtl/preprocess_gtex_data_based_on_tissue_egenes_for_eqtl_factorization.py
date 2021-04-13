import numpy as np 
import os
import sys
import pdb
import random



def extract_eqtl_factorization_tests(cross_tissue_eqtl_results_file, seed, nominal_pvalue_thresh, num_genes):
	np.random.seed(seed)
	dicti = {}
	binary_arr = []
	temp_dicti = {}
	f = open(cross_tissue_eqtl_results_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Extract relevent fields
		gene_id = data[1]
		variant_id = data[0]
		pvalue = float(data[4])
		if pvalue > nominal_pvalue_thresh:
			continue
		if gene_id not in temp_dicti:
			temp_dicti[gene_id] = []
		temp_dicti[gene_id].append(variant_id)
	f.close()

	# Randomly select which genes to use
	valid_genes = temp_dicti.keys()
	test_genes = np.random.choice(valid_genes,size=num_genes, replace=False)

	# Randomly select which variant to use for each gene
	for test_gene in test_genes:
		valid_variants = temp_dicti[test_gene]
		test_variant = np.random.choice(valid_variants, size=1, replace=False)[0]
		test_name = test_variant + ':' + test_gene
		dicti[test_name] = 1

	f = open(cross_tissue_eqtl_results_file)
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')

		gene_id = data[1]
		variant_id = data[0]
		pvalue = float(data[4])
		test_name = variant_id + ':' + gene_id
		if test_name in dicti:
			binary_arr.append(1)
			if pvalue > nominal_pvalue_thresh:
				print('assumption error')
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

# Extract array of tests to use
# Extract a dictionary for variants that maps from each variant in tests to an initialized array of length== length(sample_names)
# Extract a dictionary for genes that maps from each genes in tests to an initialized array of length== length(sample_names)
def extract_tests(tissues, gtex_egene_dir, genes_tested_in_all_tissues, valid_variants, valid_tests, num_tests):
	num_tests_per_tissue = int(round(float(num_tests)/len(tissues)))
	#num_tests_per_tissue = 1000
	# Initailize output objects
	tests = []
	variants = {}
	genes = {}
	test_dicti_final = {}
	# Loop through tissues
	for tissue in tissues:
		possible_tests = []
		# Get egene file for this tissue
		egene_file = gtex_egene_dir + tissue + '.v8.egenes.txt'
		# Stream egene file for this tissue
		f = open(egene_file)
		head_count = 0
		data = 'hi'
		for line in f:
			line = line.rstrip()
			prev_data = data
			data = line.split('\t')
			if len(data) != 33:
				continue
			# skip header
			if head_count == 0:
				head_count = head_count + 1
				header = data
				continue
			ensamble_id = data[0]
			# Skip genes not tested in all tissues
			if ensamble_id not in genes_tested_in_all_tissues:
				continue
			snp_id = data[11]
			# Skip variant if not in our list
			if snp_id not in valid_variants:
				continue
			if ensamble_id + ':' + snp_id not in valid_tests:
				continue
			# limit to autosomal chromosomes
			if snp_id.split('_')[0] == 'chrX' or snp_id.split('_')[0] == 'chrY':
				continue
			rs_id = data[18]
			qval = float(data[28])
			if qval < .05:
				possible_tests.append(ensamble_id + ':' + snp_id)
		f.close()
		if len(possible_tests) < num_tests_per_tissue:
			print('assumption error')
			pdb.set_trace()
		tests_in_this_tissue = random.sample(possible_tests, num_tests_per_tissue)
		for test in tests_in_this_tissue:
			ensamble_id = test.split(':')[0]
			snp_id = test.split(':')[1]
			tests.append(ensamble_id + ':' + snp_id)
			test_dicti_final[ensamble_id + ':' + snp_id] = 1
	return np.unique(tests), test_dicti_final

def extract_valid_variants_and_genes_and_tests(all_test_names_file):
	valid_variants = {}
	valid_genes = {}
	valid_tests = {}
	f = open(all_test_names_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		variant_id = data[1]
		valid_variants[variant_id] = 1
		valid_genes[gene_id] = 1
		valid_tests[gene_id + ':' + variant_id] = 1
	f.close()
	return valid_variants, valid_genes, valid_tests

def get_tissues(tissues_file):
	f = open(tissues_file)
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		arr.append(data[0])
	f.close()
	return arr

def extract_test_eqtl_binary_arr(all_test_names_file, test_dicti):
	binary_arr = []
	f = open(all_test_names_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		variant_id = data[1]
		test_name = gene_id + ':' + variant_id
		if test_name not in test_dicti:
			binary_arr.append(0)
		else:
			binary_arr.append(1)
	f.close()
	return np.asarray(binary_arr)

######################
# Command line args
######################
working_dir = sys.argv[1]
num_genes = int(sys.argv[2])
seed = int(sys.argv[3])
tissues_file = sys.argv[4]
gtex_egene_dir = sys.argv[5]


# Input files
# cross_tissue_eqtl_results_file = working_dir + 'cross_tissue_re_eqtl_results_merged.txt'


# Find tests to input into eqtl factorization model
# Basically we are going to limit to one snp per gene.
# and the snp-gene pairs must have pvalue < nominal_pvalue_thresh
# Randomly select which genes to use (limit to num_genes)
# and ranomly select which variants for each gene
# test_eqtl_dicti, test_eqtl_binary_arr = extract_eqtl_factorization_tests(cross_tissue_eqtl_results_file, seed, nominal_pvalue_thresh, num_genes)

tissues = get_tissues(tissues_file)

all_test_names_file = working_dir + 'all_tests.txt'
valid_variants, valid_genes, valid_tests = extract_valid_variants_and_genes_and_tests(all_test_names_file)


test_arr, test_dicti = extract_tests(tissues, gtex_egene_dir, valid_genes, valid_variants, valid_tests, num_genes)

test_eqtl_binary_arr = extract_test_eqtl_binary_arr(all_test_names_file, test_dicti)


# Generate eqtl factorization test names file
all_test_names_file = working_dir + 'all_tests.txt'
eqtl_factorization_test_names_file = working_dir + 'eqtl_factorization_tissue_egenes_input_test_names.txt'
generate_eqtl_factorization_test_names_file(all_test_names_file, eqtl_factorization_test_names_file, test_eqtl_binary_arr)

# Generate eqtl factorization gene expression file
all_gene_expression_file = working_dir + 'cross_tissue_eqtl_expression_input.txt'
eqtl_factorization_expression_file = working_dir + 'eqtl_factorization_tissue_egenes_input_expression.txt'
generate_eqtl_factorization_expression_file(all_gene_expression_file, eqtl_factorization_expression_file, test_eqtl_binary_arr)
save_as_npy_file(eqtl_factorization_expression_file)


# Generate eqtl factorization genotype expression file
all_test_genotype_file = working_dir + 'cross_tissue_eqtl_genotype_input.txt'
eqtl_factorization_genotype_file = working_dir + 'eqtl_factorization_tissue_egenes_input_genotype.txt'
generate_eqtl_factorization_genotype_file(all_test_genotype_file, eqtl_factorization_genotype_file, test_eqtl_binary_arr)
save_as_npy_file(eqtl_factorization_genotype_file)


# Add intercept to covariate file
covariate_file = working_dir + 'cross_tissue_eqtl_covariate_input.txt'
covariate_with_intercept_file = working_dir + 'cross_tissue_eqtl_covariate_w_intercept_input.txt'
add_intercept_to_covariate_file(covariate_file, covariate_with_intercept_file)