import numpy as np 
import os
import sys
import pdb




def extract_eqtl_factorization_tests(test_names_file, cross_tissue_genome_wide_sig_results_file, num_genes, ensamble_id_to_gene_symbol_id):
	dicti = {}
	binary_arr = []
	temp_dicti = {}
	f = open(cross_tissue_genome_wide_sig_results_file)
	used_genes = {}
	used_variants = {}
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
		gene_symbol_id = ensamble_id_to_gene_symbol_id[gene_id]
		if gene_symbol_id.startswith('HLA'):
			continue
		if gene_id in used_genes:
			continue
		if variant_id in used_variants:
			continue
		used_genes[gene_id] = 1
		used_variants[variant_id] = 1
		test_name = variant_id + ':' + gene_id
		if test_name in dicti:
			print('assumption error')
			pdb.set_trace()
		if len(dicti) < num_genes:
			dicti[test_name] = 1
	f.close()

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
		test_name = variant_id + ':' + gene_id
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

def generate_covariate_file(covariate_file, output_file):
	cov1 = np.loadtxt(covariate_file)
	num_samples = cov1.shape[0]
	standardized_cov = standardize_columns(cov1)
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

def generate_temp_covariate_file(output_file, pcs, geno_pcs):
	cov = np.hstack((pcs, geno_pcs))
	standardized_cov = standardize_columns(cov)
	cov_plus_intercept = np.hstack((np.ones((cov.shape[0], 1)), standardized_cov))
	np.savetxt(output_file, cov_plus_intercept, fmt="%s", delimiter='\t')

def extract_geno_pcs(geno_pc_file):
	f = open(geno_pc_file)
	head_count = 0
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			num_fields = len(data)
			continue
		if len(data) != num_fields:
			continue
		geno_pc1 = float(data[-2])
		geno_pc2 = float(data[-1])
		arr.append(np.asarray([geno_pc1, geno_pc2]))
	f.close()
	return np.asarray(arr)

def create_mapping_from_ensamble_id_to_gene_symbol_id(expr_file):
	mapping = {}
	f = open(expr_file)
	used_symbols = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		full_gene = data[0]
		ensamble_id = full_gene.split('_')[0]
		gene_symbol_id = full_gene.split('_')[1]
		if gene_symbol_id in used_symbols:
			print('assumption eororororororor')
			pdb.set_trace()
		if ensamble_id in mapping:
			print('assumption erroror')
			pdb.set_trace()
		used_symbols[gene_symbol_id] = 1
		mapping[ensamble_id] = gene_symbol_id
	f.close()
	return mapping


######################
# Command line args
######################
eqtl_input_dir = sys.argv[1]
eqtl_results_dir = sys.argv[2]
output_dir = sys.argv[3]
processed_expression_dir = sys.argv[4]

output_dir = output_dir + 'no_repeats_no_hla_'

preprocess_versions = ['scanpy_1000', 'scanpy_2000', 'scanpy_3000', 'scanpy_4000']

for preprocess_version in preprocess_versions:

	expr_file = processed_expression_dir + 'standardized_10_cap_normalized_expression_' + preprocess_version + '_hvg_all_genotyped_cells.txt'
	ensamble_id_to_gene_symbol_id = create_mapping_from_ensamble_id_to_gene_symbol_id(expr_file)

	num_genes = 100000
	# Input files
	latent_factor_interaction_eqtl_genome_wide_sig_results_file = eqtl_results_dir + 'sc_standard_eqtl_analysis_' + preprocess_version + '_hvg_eqtl_results_genome_wide_signficant_bf_fdr_0.05.txt'
	latent_factor_interaction_test_names_file = eqtl_input_dir + preprocess_version + '_hvg_eqtl_input_variant_gene_pairs.txt'

	test_eqtl_dicti, test_eqtl_binary_arr = extract_eqtl_factorization_tests(latent_factor_interaction_test_names_file, latent_factor_interaction_eqtl_genome_wide_sig_results_file, num_genes, ensamble_id_to_gene_symbol_id)

	# Generate eqtl factorization test names file
	all_test_names_file = eqtl_input_dir + preprocess_version + '_hvg_eqtl_input_variant_gene_pairs.txt'
	eqtl_factorization_test_names_file = output_dir + 'eqtl_factorization_standard_eqtl_' + preprocess_version + '_hvg_test_names.txt'
	generate_eqtl_factorization_test_names_file(all_test_names_file, eqtl_factorization_test_names_file, test_eqtl_binary_arr)

	# Generate eqtl factorization gene expression file
	all_gene_expression_file = eqtl_input_dir + preprocess_version + '_hvg_eqtl_input_expression.txt'
	eqtl_factorization_expression_file = output_dir + 'eqtl_factorization_standard_eqtl_' + preprocess_version + '_hvg_expression.txt'
	generate_eqtl_factorization_expression_file(all_gene_expression_file, eqtl_factorization_expression_file, test_eqtl_binary_arr)
	save_as_npy_file(eqtl_factorization_expression_file)


	# Generate eqtl factorization genotype expression file
	all_test_genotype_file = eqtl_input_dir + preprocess_version + '_hvg_eqtl_input_genotype.txt'
	eqtl_factorization_genotype_file = output_dir + 'eqtl_factorization_standard_eqtl_' + preprocess_version + '_hvg_unnormalized_genotype.txt'
	generate_eqtl_factorization_genotype_file(all_test_genotype_file, eqtl_factorization_genotype_file, test_eqtl_binary_arr)
	save_as_npy_file(eqtl_factorization_genotype_file)

	# Generate covariate file
	covariate_file = eqtl_input_dir + preprocess_version + '_hvg_eqtl_input_covariates.txt'
	covariate_with_intercept_file = output_dir + 'eqtl_factorization_standard_eqtl_' + preprocess_version + '_hvg_covariates.txt'
	generate_covariate_file(covariate_file, covariate_with_intercept_file)

	# Generate sample overlap file
	sample_overlap_input_file = eqtl_input_dir + preprocess_version + '_hvg_eqtl_input_sample_overlap.txt'
	sample_overlap_output_file = output_dir + 'eqtl_factorization_standard_eqtl_' + preprocess_version + '_hvg_sample_overlap.txt'
	generate_sample_overlap_file(sample_overlap_input_file, sample_overlap_output_file)



