import numpy as np 
import os
import sys
import pdb




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


######################
# Command line args
######################
working_dir = sys.argv[1]
num_genes = int(sys.argv[2])
seed = int(sys.argv[3])
nominal_pvalue_thresh = float(sys.argv[4])


# Input files
cross_tissue_eqtl_results_file = working_dir + 'cross_tissue_re_eqtl_results_merged.txt'


# Find tests to input into eqtl factorization model
# Basically we are going to limit to one snp per gene.
# and the snp-gene pairs must have pvalue < nominal_pvalue_thresh
# Randomly select which genes to use (limit to num_genes)
# and ranomly select which variants for each gene
test_eqtl_dicti, test_eqtl_binary_arr = extract_eqtl_factorization_tests(cross_tissue_eqtl_results_file, seed, nominal_pvalue_thresh, num_genes)


# Generate eqtl factorization test names file
all_test_names_file = working_dir + 'all_tests.txt'
eqtl_factorization_test_names_file = working_dir + 'eqtl_factorization_p_thresh_' + str(nominal_pvalue_thresh) + '_input_test_names.txt'
generate_eqtl_factorization_test_names_file(all_test_names_file, eqtl_factorization_test_names_file, test_eqtl_binary_arr)

# Generate eqtl factorization gene expression file
all_gene_expression_file = working_dir + 'cross_tissue_eqtl_expression_input.txt'
eqtl_factorization_expression_file = working_dir + 'eqtl_factorization_p_thresh_' + str(nominal_pvalue_thresh) + '_input_expression.txt'
generate_eqtl_factorization_expression_file(all_gene_expression_file, eqtl_factorization_expression_file, test_eqtl_binary_arr)
save_as_npy_file(eqtl_factorization_expression_file)


# Generate eqtl factorization genotype expression file
all_test_genotype_file = working_dir + 'cross_tissue_eqtl_genotype_input.txt'
eqtl_factorization_genotype_file = working_dir + 'eqtl_factorization_' + str(nominal_pvalue_thresh) + '_input_genotype.txt'
generate_eqtl_factorization_genotype_file(all_test_genotype_file, eqtl_factorization_genotype_file, test_eqtl_binary_arr)
save_as_npy_file(eqtl_factorization_genotype_file)



