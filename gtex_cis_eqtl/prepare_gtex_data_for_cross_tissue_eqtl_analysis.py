import numpy as np 
import os
import sys
import pdb






def order_variant_gene_pairs_file_by_chromosome(all_unordered_variant_gene_pairs_file, all_variant_gene_pairs_file):
	t = open(all_variant_gene_pairs_file, 'w')
	t.write('gene_id\tvariant_id\n')
	for chrom_num in range(1,23):
		chrom_string = 'chr' + str(chrom_num)
		f = open(all_unordered_variant_gene_pairs_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			variant_id = data[1]
			chrom_variant_on = variant_id.split('_')[0]
			if chrom_variant_on != chrom_string:
				continue
			t.write(line + '\n')
		f.close()
	t.close()


def extract_ordered_samples_according_to_covariate_file(covariate_file):
	data = np.loadtxt(covariate_file, dtype=str, delimiter='\t')
	return data[1:,0]

def extract_genes_on_chromosome_for_testing(all_variant_gene_pairs_file, chrom_num):
	chrom_string = 'chr' + str(chrom_num)
	f = open(all_variant_gene_pairs_file)
	arr = []
	dicti = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_id = data[1]
		chrom_variant_on = data[1].split('_')[0]
		if chrom_variant_on != chrom_string:
			continue
		gene_id = data[0]
		arr.append(gene_id)
		dicti[gene_id] = np.zeros(1)
	f.close()
	return arr, dicti

def extract_variants_on_chromosome_for_testing(all_variant_gene_pairs_file, chrom_num, ordered_samples):
	chrom_string = 'chr' + str(chrom_num)
	f = open(all_variant_gene_pairs_file)
	arr = []
	dicti = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_id = data[1]
		chrom_variant_on = data[1].split('_')[0]
		if chrom_variant_on != chrom_string:
			continue
		gene_id = data[0]
		arr.append(variant_id)
		dicti[variant_id] = np.zeros(len(ordered_samples))
	f.close()
	return arr, dicti

def fill_in_mapping_from_gene_name_to_expression_vector(genes_dicti, expression_file, ordered_samples):
	f = open(expression_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Header
		if head_count == 0:
			head_count = head_count + 1
			# Simple error checking to ensure samples in expression file are in correct order
			if np.array_equal(data[1:], ordered_samples) == False:
				print('assumption eroror')
				pdb.set_trace()
			continue
		# standard line
		gene_id = data[0]
		if gene_id not in genes_dicti:
			continue
		# Simple error checking
		if len(genes_dicti[gene_id]) != 1:
			print('assumption eorror')
			pdb.set_trace()
		genes_dicti[gene_id] = data[1:]
	f.close()
	return genes_dicti


def generate_eqtl_expression_input_data(all_variant_gene_pairs_file, ordered_samples, expression_file, eqtl_expression_file):
	# Open output file handle
	t = open(eqtl_expression_file, 'w')
	# Do this one chromosome at a time (for memory purposes)
	for chrom_num in range(1, 23):
		print(chrom_num)
		# Extract ordered list of genes of this chromosome and dictionary list of genes on this chromosomes
		ordered_genes_arr, genes_dicti = extract_genes_on_chromosome_for_testing(all_variant_gene_pairs_file, chrom_num)

		# Fill up genes_dicti with expression data
		genes_dicti = fill_in_mapping_from_gene_name_to_expression_vector(genes_dicti, expression_file, ordered_samples)

		for gene_id in ordered_genes_arr:
			# Simple error checking
			if len(genes_dicti[gene_id]) != len(ordered_samples):
				print('assumption eoror')
				pdb.set_trace()
			t.write('\t'.join(genes_dicti[gene_id]) + '\n')
	t.close()


# Fill in 'variants' dictionary with genotype values
def add_genotype_values_to_data_structure(variants, sample_names, gtex_genotype_dir, chrom_num):
	used_variants = {}
	genotype_file = gtex_genotype_dir + 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + str(chrom_num) + '_dosage_MAF_05.txt'
	head_count = 0
	# Stream genotype file for this chromosome
	f = open(genotype_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			# Header contains indi_ids
			indi_ids = data
			# Create mapping from column index to sample inidices
			mapping = {}
			for index, indi_id in enumerate(indi_ids):
				valid = False
				arr = []
				for index2, sample_name in enumerate(sample_names):
					if sample_name.split(':')[0] == indi_id:
						valid = True
						arr.append(index2)
				if valid == True:
					mapping[index] = np.asarray(arr)
			continue
		snp_id = data[0]
		miss_counter = 0
		# skip lines where variant not in 'variants' dictionary
		if snp_id not in variants:
			continue
		try:
			genotype_vec = np.asarray(data[1:]).astype(float)
		except:
			miss_counter = miss_counter + 1
			for index, genotype in enumerate(data[1:]):
				if index not in mapping or genotype == '-':
					continue
				variants[snp_id][mapping[index]] = float(genotype)
			continue
		used_variants[snp_id] = 1
		for index, genotype in enumerate(genotype_vec):
			if index not in mapping:
				continue
			variants[snp_id][mapping[index]] = genotype
	f.close()
	# Some quick error chekcing
	if np.array_equal(sorted(used_variants.keys()), sorted(variants.keys())) == False:
		print('assumption erorr')
		pdb.set_trace()
	print(str(miss_counter) + ' variants with non-float data on chrom ' + str(chrom_num))
	return variants

def generate_eqtl_genotype_input_data(all_variant_gene_pairs_file, ordered_samples, gtex_genotype_dir, eqtl_genotype_file):
	# Open output file handle
	t = open(eqtl_genotype_file, 'w')
	# Do this one chromosome at a time (for memory purposes)
	for chrom_num in range(1, 23):
		# Extract ordered list of genes of this chromosome and dictionary list of genes on this chromosomes
		ordered_variant_arr, variant_dicti = extract_variants_on_chromosome_for_testing(all_variant_gene_pairs_file, chrom_num, ordered_samples)
		# Fill in 'variants' dictionary with genotype values
		variant_dicti = add_genotype_values_to_data_structure(variant_dicti, ordered_samples, gtex_genotype_dir, chrom_num)
		# Print to output file
		for variant_id in ordered_variant_arr:
			# Simple error checking
			if len(variant_dicti[variant_id]) != len(ordered_samples):
				print('assumption eoror')
				pdb.set_trace()
			# print to output
			t.write('\t'.join(variant_dicti[variant_id].astype(str)) + '\n')
	t.close()

def generate_eqtl_covariate_input_data(covariate_file, eqtl_covariate_file):
	data = np.loadtxt(covariate_file, dtype=str, delimiter='\t')
	np.savetxt(eqtl_covariate_file, data[1:,1:], fmt="%s", delimiter='\t')

def generate_eqtl_interaction_factor_input_data(covariate_file, eqtl_interaction_factor_file):
	data = np.loadtxt(covariate_file, dtype=str, delimiter='\t')
	header = data[0,:]
	# get expression pc column indices
	indices = []
	for i, ele in enumerate(header):
		if ele.startswith('PC'):
			indices.append(i)
	indices = np.asarray(indices)
	# Save to output file
	np.savetxt(eqtl_interaction_factor_file, data[1:, indices], fmt="%s", delimiter='\t')


def extract_categorical_feature_matrix(raw_feature_cov):
	all_features = []
	for raw_feature_str in raw_feature_cov:
		all_features.append(raw_feature_str)
	unique_categories = np.unique(np.asarray(all_features))
	mapping = {}
	for i, category in enumerate(unique_categories):
		mapping[category] = i

	cat_feature = np.zeros((len(raw_feature_cov), len(unique_categories) -1))

	for i, raw_feature_string in enumerate(raw_feature_cov):
		column_num = mapping[raw_feature_string]
		if column_num == (len(unique_categories) -1):
			continue
		cat_feature[i, column_num] = 1.0
	return(cat_feature)

def generate_eqtl_covariate_includinging_tissue_identity_input_data(residual_expression_covariate_file, eqtl_residual_covariate_file):
	data = np.loadtxt(residual_expression_covariate_file, dtype=str, delimiter='\t')
	standard_cov = data[1:,1:]
	# Extract sample tissue names
	sample_ids = data[1:,0]
	sample_tissue_names = []
	for sample_id in sample_ids:
		sample_tissue_names.append(sample_id.split(':')[1])
	sample_tissue_names = np.asarray(sample_tissue_names)
	# Convert sample tissue names vector to one-hot encoding binary matrix
	tissue_mat = extract_categorical_feature_matrix(sample_tissue_names)
	# create concatenated covariate matrix
	full_cov = np.hstack((standard_cov, tissue_mat.astype(str)))
	np.savetxt(eqtl_residual_covariate_file, full_cov, fmt="%s", delimiter='\t')

def generate_eqtl_interaction_factors_from_tissue_identity(residual_expression_covariate_file, eqtl_known_tissue_interaction_file):
	data = np.loadtxt(residual_expression_covariate_file, dtype=str, delimiter='\t')
	# Extract sample tissue names
	sample_ids = data[1:,0]
	sample_tissue_names = []
	for sample_id in sample_ids:
		sample_tissue_names.append(sample_id.split(':')[1])
	sample_tissue_names = np.asarray(sample_tissue_names)
	# Convert sample tissue names vector to one-hot encoding binary matrix
	tissue_mat = extract_categorical_feature_matrix(sample_tissue_names)
	np.savetxt(eqtl_known_tissue_interaction_file, tissue_mat.astype(str), fmt="%s", delimiter='\t')

#####################
# Command Line args
######################
tissues_file = sys.argv[1]
gtex_genotype_dir = sys.argv[2]
output_dir = sys.argv[3]




# Already generated data
all_unordered_variant_gene_pairs_file = output_dir + 'all_tests_unordered.txt'
expression_file = output_dir + 'cross_tissue_tpm_standardized.txt'
covariate_file = output_dir + 'covariates.txt'


# Re-order variant gene pairs file
all_variant_gene_pairs_file = output_dir + 'all_tests.txt'
order_variant_gene_pairs_file_by_chromosome(all_unordered_variant_gene_pairs_file, all_variant_gene_pairs_file)



# Extract array of ordered samples (according to already generated covariate file)
ordered_samples = extract_ordered_samples_according_to_covariate_file(covariate_file)

# Save matrix of covariates of dimension num_samplesXnum_covariates
eqtl_covariate_file = output_dir + 'cross_tissue_eqtl_covariate_input.txt'
generate_eqtl_covariate_input_data(covariate_file, eqtl_covariate_file)


# Extract expression data for eqtl analysis
eqtl_expression_file = output_dir + 'cross_tissue_eqtl_expression_input.txt'
generate_eqtl_expression_input_data(all_variant_gene_pairs_file, ordered_samples, expression_file, eqtl_expression_file)


# Extract genotype data for eqtl analysis
eqtl_genotype_file = output_dir + 'cross_tissue_eqtl_genotype_input.txt'
generate_eqtl_genotype_input_data(all_variant_gene_pairs_file, ordered_samples, gtex_genotype_dir, eqtl_genotype_file)



