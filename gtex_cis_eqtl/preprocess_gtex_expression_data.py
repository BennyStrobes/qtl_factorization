from sklearn import linear_model
import numpy as np 
import os
import sys
import pdb
import gzip
import random
import pandas as pd
import rnaseqnorm
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

def get_tissues(file_name):
	f = open(file_name)
	arr = []
	arr2 = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		arr.append(data[0])
		arr2.append(data[1])
	f.close()
	return arr, arr2

def regress_out_covariates_and_center(genotype_input_file, covariate_file, genotype_output_file):
	# Load in covariates
	covariate_raw = np.transpose(np.loadtxt(covariate_file, dtype=str, delimiter='\t'))
	covariate_names = covariate_raw[1:,0]
	covariate_samples = covariate_raw[0,1:]
	covariate_mat = covariate_raw[1:,1:].astype(float)
	
	# Load in expression data
	genotype_raw = np.loadtxt(genotype_input_file, dtype=float, delimiter='\t')
	genotype_mat = np.transpose(genotype_raw)

	# Initialize output matrix 
	num_samples = genotype_mat.shape[0]
	num_genes = genotype_mat.shape[1]

	t = open(genotype_output_file, 'w')


	model = linear_model.LinearRegression(fit_intercept=True) 
	modelfit = model.fit(np.transpose(covariate_mat),genotype_mat)
	pred = modelfit.predict(np.transpose(covariate_mat))

	resid = genotype_mat - pred
	for gene_number in range(num_genes):
		# print(np.std(resid[:,gene_number]))
		#residual_expression = regress_out_covariates_for_one_gene(expr_mat[:,gene_number], covariate_mat)
		#gene_id = gene_ids[gene_number]
		t.write('\t'.join(resid[:,gene_number].astype(str)) + '\n')
	t.close()

def regress_out_covariates_and_standardize(genotype_input_file, covariate_file, genotype_output_file):
	# Load in covariates
	covariate_raw = np.transpose(np.loadtxt(covariate_file, dtype=str, delimiter='\t'))
	covariate_names = covariate_raw[1:,0]
	covariate_samples = covariate_raw[0,1:]
	covariate_mat = covariate_raw[1:,1:].astype(float)
	
	# Load in expression data
	genotype_raw = np.loadtxt(genotype_input_file, dtype=float, delimiter='\t')
	genotype_mat = np.transpose(genotype_raw)


	# Initialize output matrix 
	num_samples = genotype_mat.shape[0]
	num_genes = genotype_mat.shape[1]


	for gene_num in range(num_genes):
		genotype_mat[:, gene_num] = (genotype_mat[:, gene_num] - np.mean(genotype_mat[:, gene_num]))/np.std(genotype_mat[:, gene_num])

	t = open(genotype_output_file, 'w')


	model = linear_model.LinearRegression(fit_intercept=True) 
	modelfit = model.fit(np.transpose(covariate_mat),genotype_mat)
	pred = modelfit.predict(np.transpose(covariate_mat))

	resid = genotype_mat - pred
	for gene_number in range(num_genes):
		# print(np.std(resid[:,gene_number]))
		#residual_expression = regress_out_covariates_for_one_gene(expr_mat[:,gene_number], covariate_mat)
		#gene_id = gene_ids[gene_number]
		t.write('\t'.join(resid[:,gene_number].astype(str)) + '\n')
	t.close()





def regress_out_covariates(expression_input_file, covariate_file, expression_output_file):
	# Load in covariates
	covariate_raw = np.transpose(np.loadtxt(covariate_file, dtype=str, delimiter='\t'))
	covariate_names = covariate_raw[1:,0]
	covariate_samples = covariate_raw[0,1:]
	covariate_mat = covariate_raw[1:,1:].astype(float)
	
	# Load in expression data
	expression_raw = np.loadtxt(expression_input_file, dtype=str, delimiter='\t')
	expression_samples = expression_raw[1:,0]
	gene_ids = expression_raw[0,1:]
	expr_mat = expression_raw[1:,1:].astype(float)
	# Simple error checking
	if np.array_equal(expression_samples, covariate_samples) == False:
		print('assumption error!')
		pdb.set_trace()

	# Initialize output matrix 
	num_samples = expr_mat.shape[0]
	num_genes = expr_mat.shape[1]

	t = open(expression_output_file, 'w')
	t.write('Gene_id\t' + '\t'.join(expression_samples) + '\n')


	model = linear_model.LinearRegression(fit_intercept=True) 
	modelfit = model.fit(np.transpose(covariate_mat),expr_mat)
	pred = modelfit.predict(np.transpose(covariate_mat))

	resid = expr_mat - pred
	for gene_number in range(num_genes):
		# print(np.std(resid[:,gene_number]))
		#residual_expression = regress_out_covariates_for_one_gene(expr_mat[:,gene_number], covariate_mat)
		gene_id = gene_ids[gene_number]
		t.write(gene_id + '\t' + '\t'.join(resid[:,gene_number].astype(str)) + '\n')
	t.close()
	'''
	# Open expression input and output files
	f = gzip.open(expression_input_file)
	t = open(expression_output_file, 'w')
	# Stream input file line by line
	head_count = 0  # to identify header line
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			# check to make sure expression samples are in the same order as covariate sampels
			if np.array_equal(data[4:],covariate_samples) == False:
				print('Assumption error: samples not alligned')
			# Print header to output file
			t.write('\t'.join(data[3:]) + '\n')
			continue
		gene_id = data[3]
		expression = np.asarray(data[4:]).astype(float)
		residual_expression = regress_out_covariates_for_one_gene(expression, covariate_mat)
		t.write(gene_id + '\t' + '\t'.join(residual_expression.astype(str)) + '\n')
	f.close()
	t.close()
	'''

def get_sample_names(tissues, gtex_expression_dir, sample_name_file, gtex_individual_information_file):
	# Initialize arr
	samples = []
	# get samples in tisssue
	for tissue in tissues:
		expression_file = gtex_expression_dir + tissue + '.v8.normalized_expression.bed.gz'
		f = gzip.open(expression_file)
		# sample names are only in the header
		head_count = 0
		for line in f:
			if head_count == 0:
				line = line.rstrip()
				data = line.split()
				head_count = head_count +1
				for indi_id in data[4:]:
					samples.append(indi_id + ':' + tissue)
				continue
			break
		f.close()
	indis = {}
	for sample in samples:
		indi = sample.split(':')[0]
		indis[indi] = 1
	# Get mapping from sample_name to index
	sample_to_index = {}
	# Print to output file
	t = open(sample_name_file,'w')
	for i,sample in enumerate(samples):
		t.write(sample + '\n')
		sample_to_index[sample] = i
	t.close()
	# Extract covariates for these samples
	dicti = {}
	f = open(gtex_individual_information_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# header
		if head_count == 0:
			header = np.copy(data)
			head_count = head_count +1
			continue
		indi_id = data[0]
		cohort = data[1]
		sex = data[2]
		age = data[3]
		race = data[4]
		ethnicity = data[5]
		height = data[6]
		weight = data[8]
		bmi = data[10]
		# Extract ischemic time in minutes
		ischemic_time_string = data[12]
		ischemic_time_data = ischemic_time_string.split(',')
		ischemic_hours = float(ischemic_time_string.split(' hour')[0])
		ischemic_minutes = float(ischemic_time_string.split(',')[1].split(' minute')[0])
		total_ischemic_minutes = ischemic_hours*60.0 + ischemic_minutes		
		# error checking
		if data[7] != 'in' or data[9] != 'lb':
			print('assumptino erorror')
			pdb.set_trace()
		dicti[indi_id] = [cohort, sex, age, race, height, weight, bmi, str(total_ischemic_minutes)]
	f.close()
	# Some quick error checking
	for indi_id in dicti.keys():
		if len(dicti[indi_id]) != 8:
			print('assumption erorr')
			pdb.set_trace()
	return samples, sample_to_index, dicti

# Extract array of tests to use
# Extract a dictionary for variants that maps from each variant in tests to an initialized array of length== length(sample_names)
# Extract a dictionary for genes that maps from each genes in tests to an initialized array of length== length(sample_names)
def extract_tests(tissues, gtex_egene_dir, num_samples, genes_tested_in_all_tissues, valid_variants):
	#num_tests_per_tissue = int(round(25000.0/len(tissues)))
	num_tests_per_tissue = 1000
	# Initailize output objects
	tests = []
	variants = {}
	genes = {}
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
			variants[snp_id] = np.zeros(num_samples)
			genes[ensamble_id] = np.zeros(num_samples)
	return np.unique(tests), variants, genes

# Print test names to output file
def print_test_names(tests, test_name_file):
	f = open(test_name_file, 'w')
	for test in tests:
		f.write(test + '\n')
	f.close()

# We are going to limit analysis to genes tested in all tissues
def get_genes_tested_in_all_tissues(tissues, gtex_expression_dir):
	# Initialize dictionaries to keep track of which genes were in all tissues
	gene_counts = {}
	valid_genes = {}
	# For each tissue, keep track fo which genes were used
	for tissue in tissues:
		gene_file_name = gtex_expression_dir + tissue + '.v8.normalized_expression.bed.gz'
		head_count = 0
		f = gzip.open(gene_file_name)
		for line in f:
			line = line.rstrip()
			data = line.split()
			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue
			ensamble_id = data[3]
			# Add gene to dictionary if not observed
			if ensamble_id not in gene_counts:
				gene_counts[ensamble_id] = 0
			# Keep track of how many tissues this gene was observed in 
			gene_counts[ensamble_id] = gene_counts[ensamble_id] + 1
		f.close()
	# Loop through all observed ensamble ids
	# Only take genes that are expressed in all tissues
	for ensamble_id in gene_counts.keys():
		if gene_counts[ensamble_id] == len(tissues):
			valid_genes[ensamble_id] = 1
	return valid_genes

# Fill in 'genes' dictionary with expression values from each tissue
def add_expression_values_to_data_structure(expr_file, sample_to_index, genes):
	f = open(expr_file)
	# to identify header
	head_count = 0
	# Stream file
	for line in f:
		line = line.rstrip()
		data = line.split()
		# header
		if head_count == 0:
			head_count = head_count + 1
			# Creating mapping from index to sample name
			mapping = {}
			for i, indi_id in enumerate(data[1:]):
				mapping[i] = indi_id
			continue
		ensamble_id = data[0]
		# Only consider lines where ensamble_id is in genes dictionary
		if ensamble_id not in genes:
			continue
		expr_vec = np.asarray(data[1:]).astype(float)
		for index, expr in enumerate(expr_vec):
			sample_name = mapping[index] 
			genes[ensamble_id][sample_to_index[sample_name]] = expr
	f.close()
	return genes

# Fill in 'genes' dictionary with expression values from each tissue
def add_expression_values_to_data_structure_t(expr_file, sample_to_index, genes):
	aa = np.loadtxt(expr_file,dtype=str,delimiter='\t')
	np.savetxt('temp.txt', np.transpose(aa), fmt="%s", delimiter='\t')
	f = open('temp.txt')
	# to identify header
	head_count = 0
	# Stream file
	for line in f:
		line = line.rstrip()
		data = line.split()
		# header
		if head_count == 0:
			head_count = head_count + 1
			# Creating mapping from index to sample name
			mapping = {}
			for i, indi_id in enumerate(data[1:]):
				mapping[i] = indi_id
			continue
		ensamble_id = data[0]
		# Only consider lines where ensamble_id is in genes dictionary
		if ensamble_id not in genes:
			continue
		expr_vec = np.asarray(data[1:]).astype(float)
		for index, expr in enumerate(expr_vec):
			sample_name = mapping[index] 
			genes[ensamble_id][sample_to_index[sample_name]] = expr
	f.close()
	return genes


# Fill in 'variants' dictionary with genotype values
def add_genotype_values_to_data_structure(variants, sample_names, gtex_genotype_dir):
	used_variants = {}
	# loop through chromosomes
	for chrom_num in range(1,23):
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
			# skip lines where variant not in 'variants' dictionary
			if snp_id not in variants:
				continue
			try:
				genotype_vec = np.asarray(data[1:]).astype(float)
			except:
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
	return variants


# Fill in 'variants' dictionary with genotype values
def get_variants_we_have_genotype_data_for( gtex_genotype_dir):
	used_variants = {}
	# loop through chromosomes
	for chrom_num in range(1,23):
		genotype_file = gtex_genotype_dir + 'GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_chr' + str(chrom_num) + '_dosage_MAF_05.txt'
		head_count = 0
		# Stream genotype file for this chromosome
		f = open(genotype_file)
		for line in f:
			line = line.rstrip()
			data = line.split()
			if head_count == 0:
				head_count = head_count + 1
				continue
			snp_id = data[0]
			valid = True
			for ele in data[1:]:
				if ele == '-':
					valid=False
			if valid == True:
				used_variants[snp_id] = data[0]
		f.close()
	return used_variants

def print_big_matrix(output_file, tests, data_struct, index):
	t = open(output_file,'w')
	for test in tests:
		name = test.split(':')[index]
		t.write('\t'.join(data_struct[name].astype(str)) + '\n')
	t.close()

def print_individual_id_file_and_tissue_file(sample_names, output_file, output_tissue_file):
	# Create mapping for indi_id to index
	counter = 0
	mapping = {}
	for sample_name in sample_names:
		indi_id = sample_name.split(':')[0]
		if indi_id not in mapping:
			mapping[indi_id] = counter
			counter = counter  + 1
	t = open(output_file, 'w')
	for sample_name in sample_names:
		indi_id = sample_name.split(':')[0]
		index = mapping[indi_id]
		t.write(str(index) + '\n')
	t.close()
	t = open(output_tissue_file, 'w')
	for sample_name in sample_names:
		tissue_name = sample_name.split(':')[1]
		t.write(tissue_name + '\n')
	t.close()


def remove_genes_with_zero_expression(tpm_matrix, ordered_genes):
	valid_indices = []
	for gene_number in range(len(ordered_genes)):
		if np.sum(tpm_matrix[:,gene_number]) != 0:
			valid_indices.append(gene_number)
	valid_indices = np.asarray(valid_indices)
	return tpm_matrix[:,valid_indices], np.asarray(ordered_genes)[valid_indices]

# Generate TPM expression matrix
def generate_tpm_expression_matrix(tissues, tissues_alt, ordered_sample_names, genes_tested_in_all_tissues, gtex_tpm_dir, output_file):
	# Generate mapping from sample_name to row number
	sample_name_to_row_number = {}
	for index, sample_name in enumerate(ordered_sample_names):
		sample_name_to_row_number[sample_name] = index
	# Generate mapping from gene to column number
	ordered_genes = sorted(genes_tested_in_all_tissues.keys())
	gene_to_column_number = {}
	for index, gene in enumerate(ordered_genes):
		gene_to_column_number[gene] = index
	# Num samples and number of genes
	num_samples = len(sample_name_to_row_number)
	num_genes = len(ordered_genes)
	# Initialize tpm matrix
	tpm_matrix = np.zeros((num_samples, num_genes))
	counter = 0
	used = {}
	# Loop through each tissue and fill in the tpm matrix
	for tissue_index, tissue in enumerate(tissues):
		tissue_alt = tissues_alt[tissue_index]
		# Stream tpm file in this tissue
		tpm_file = gtex_tpm_dir + tissue_alt + '.txt'
		head_count = 0
		f = open(tpm_file)
		for line in f:
			line = line.rstrip()
			data = line.split()
			# Header
			if head_count == 0:
				head_count = head_count + 1
				# Get ordered list of sample names
				sample_names = []
				for sample_name_temp in data[2:]:
					sample_name = sample_name_temp.split('-')[0] + '-' + sample_name_temp.split('-')[1] + ':' + tissue
					sample_names.append(sample_name)
				continue
			# Get relevent fields from line
			ensamble_id = data[0]
			# Skip genes we are not interested in
			if ensamble_id not in gene_to_column_number:
				continue
			tpm_counts = np.asarray(data[2:]).astype(float)
			# Loop through samples and add sample/tpm count to tpm_matrix
			for index, sample_name in enumerate(sample_names):
				# Ignore samples that aren't in our list
				if sample_name not in sample_name_to_row_number:
					continue
				# Get row corresponding to this smample
				row_index = sample_name_to_row_number[sample_name]
				# Get column corresponding to this gene
				column_index = gene_to_column_number[ensamble_id]

				ele_name = str(row_index) + '_' + str(column_index)

				# Add tpm count to matrix
				tpm_matrix[row_index, column_index] = tpm_counts[index]
				counter = counter + 1
		f.close()
	# Remove genes with zero expression across all samples
	filtered_tpm_matrix, filtered_orderd_genes = remove_genes_with_zero_expression(tpm_matrix, ordered_genes)
	log_filtered_tpm_matrix = np.log2(filtered_tpm_matrix + 1.0)
	# Print to output file
	t = open(output_file, 'w')
	# print header
	t.write('GeneId\t' + '\t'.join(filtered_orderd_genes) + '\n')
	for sample_num, sample_name in enumerate(ordered_sample_names):
		tpm_expr = log_filtered_tpm_matrix[sample_num,:].astype(str)
		t.write(sample_name + '\t' + '\t'.join(tpm_expr) + '\n')
	t.close()


def standardize_expression(tpm_expression_matrix_file, standardized_tpm_expression_matrix_file):
	tpm_full = np.loadtxt(tpm_expression_matrix_file, dtype=str,delimiter='\t')
	tpm = tpm_full[1:,1:].astype(float)
	samples = tpm_full[1:,0]
	genes = tpm_full[0,1:]
	# Quantile normalize the samples
	df = pd.DataFrame(np.transpose(tpm))
	#rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
	#temp_out = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
	#tpm_quantile_normalized = np.transpose(np.asarray(temp_out))
	temp_out = rnaseqnorm.normalize_quantiles(df)
	norm_df = rnaseqnorm.inverse_normal_transform(temp_out)
	standardized_tpm = np.transpose(np.asarray(norm_df))
	###
	#tpm_quantile_normalized = np.transpose(np.asarray(temp_out))
	###

	# Standardize the genes
	#num_genes = tpm_quantile_normalized.shape[1]
	#num_samples = tpm_quantile_normalized.shape[0]

	####
	#standardized_tpm = np.zeros((num_samples, num_genes))
	#for gene_num in range(num_genes):
	#	standardized_tpm[:,gene_num] = (tpm_quantile_normalized[:, gene_num] - np.mean(tpm_quantile_normalized[:, gene_num]))/np.std(tpm_quantile_normalized[:, gene_num])
	####
	# Print to output file
	t = open(standardized_tpm_expression_matrix_file, 'w')
	# print header
	t.write('SampleId\t' + '\t'.join(samples) + '\n')
	for gene_num, gene_name in enumerate(genes):
		#expr = tpm_quantile_normalized[sample_num, :].astype(str)
		###
		expr = standardized_tpm[:, gene_num].astype(str)
		###
		t.write(gene_name + '\t' + '\t'.join(expr) + '\n')
	t.close()
	'''
	t = open(standardized_tpm_expression_matrix_file, 'w')
	# print header
	t.write('GeneId\t' + '\t'.join(genes) + '\n')
	for sample_num, sample_name in enumerate(samples):
		#expr = tpm_quantile_normalized[sample_num, :].astype(str)
		###
		expr = standardized_tpm[sample_num, :].astype(str)
		###
		t.write(sample_name + '\t' + '\t'.join(expr) + '\n')
	t.close()
	'''

# First extract matrix of samples X genotype PCs
def extract_genotype_pcs(gtex_covariate_dir, sample_names, tissues, num_genotype_pcs):
	# Create dictionary of sample_name to genotype_PCs
	dicti = {}
	for sample_name in sample_names:
		dicti[sample_name] = 0
	for tissue in tissues:
		covariate_file = gtex_covariate_dir + tissue + '.v8.covariates.txt'
		cov_data = np.loadtxt(covariate_file, dtype=str,delimiter='\t')
		covs = cov_data[1:,1:].astype(float)
		cov_names = cov_data[1:,0]
		ids = cov_data[0,1:]
		for index, individual in enumerate(ids):
			individual_name = individual + ':' + tissue
			if individual_name in dicti:
				dicti[individual_name] = covs[:5,index]
	# Check
	new_dicti = {}
	for sample_name in dicti.keys():
		indi = sample_name.split(':')[0]
		pc_vec = dicti[sample_name]
		if indi not in new_dicti:
			new_dicti[indi] = pc_vec
		else:
			if np.array_equal(pc_vec, new_dicti[indi]) == False:
				print('assumption eror')
				pdb.set_trace()
	genotype_pc_mat = np.zeros((len(sample_names), num_genotype_pcs))
	for index, sample_name in enumerate(sample_names):
		genotype_pc_mat[index,:] = dicti[sample_name]
	return genotype_pc_mat


def extract_covariates(expr_file, output_file, num_expression_pcs, gtex_covariate_dir, sample_names, tissues):
	# First extract matrix of samples X genotype PCs
	num_genotype_pcs = 5
	genotype_pcs = extract_genotype_pcs(gtex_covariate_dir, sample_names, tissues, num_genotype_pcs)
	# Load in expression data
	expr_full = np.transpose(np.loadtxt(expr_file, dtype=str,delimiter='\t'))
	expr = expr_full[1:,1:].astype(float)
	samples = expr_full[1:,0]
	genes = expr_full[0,1:]

	# Compute gene expression pcs on expression data
	uuu, sss, vh = np.linalg.svd(np.transpose(expr))
	expr_pc_loadings = np.transpose(vh)[:,:num_expression_pcs]
	# print to output file 
	pc_names = []
	for pc_num in range(num_expression_pcs):
		pc_names.append('PC' + str(pc_num))
	for genotype_pc in range(num_genotype_pcs):
		pc_names.append('genotype_PC' + str(genotype_pc))
	t = open(output_file, 'w')
	t.write('Sample_name\t' + '\t'.join(pc_names) + '\n')
	for sample_num, sample_name in enumerate(samples):
		t.write(sample_name + '\t' + '\t'.join(expr_pc_loadings[sample_num,:].astype(str)))
		t.write('\t' + '\t'.join(genotype_pcs[sample_num,:].astype(str)) + '\n')
	t.close()


def get_genes_we_have_expression_data_for(file_name):
	genes = {}
	head_count = 0
	f = open(file_name)
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		genes[data[0]] = 1
	return genes


def print_sample_covariates(sample_names, individual_covariates, xcell_covariates, xcell_covariate_names, output_file):
	t = open(output_file, 'w')
	t.write('sample_id\ttissue_type\tcohort\tsex\tage\trace\theight\tweight\tbmi\tischemic_time\t' + '\t'.join(xcell_covariate_names) + '\n')
	for sample_name in sample_names:
		individual_id = sample_name.split(':')[0]
		tissue_type = sample_name.split(':')[1]
		t.write(sample_name + '\t' + tissue_type + '\t' + '\t'.join(individual_covariates[individual_id]) + '\t' + '\t'.join(xcell_covariates[sample_name]) + '\n')
	t.close()


def print_test_effect_sizes(tissues, gtex_eqtl_dir, test_name_file, test_effect_size_file):
	# Extract test names
	test_names = []
	test_names_dicti = {}
	f = open(test_name_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		test_names.append(line)
		test_names_dicti[line] = np.zeros(len(tissues)) - 400
	f.close()
	# For each tissue, fill in observed effect sizes
	for tissue_index, tissue in enumerate(tissues):
		print(tissue)
		eqtl_file = gtex_eqtl_dir + tissue + '.allpairs.txt'
		f = open(eqtl_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split()
			if head_count == 0:
				head_count = head_count + 1
				continue
			ensamble_id = data[0]
			variant_id = data[1]
			line_test_name = ensamble_id + ':' + variant_id
			if line_test_name not in test_names_dicti:
				continue
			if test_names_dicti[line_test_name][tissue_index] != -400:
				print('assumption error!')
			test_names_dicti[line_test_name][tissue_index] = float(data[7])
		f.close()
	t = open(test_effect_size_file, 'w')
	t.write('\t'.join(tissues) + '\n')
	for test_name in test_names:
		effect_size_vec = test_names_dicti[test_name]
		for index in range(len(tissues)):
			if effect_size_vec[index] == -400:
				print('assumption error!')
				pdb.set_trace()
		effect_size_vec = effect_size_vec.astype(str)
		t.write('\t'.join(effect_size_vec) + '\n')
	t.close()

def extract_all_variant_gene_pairs(tissues, gtex_eqtl_dir, tss_distance_thresh, valid_genes, valid_variants, all_tests_file):
	tests = {}
	for tissue in tissues:
		print(tissue)
		file_name = gtex_eqtl_dir + tissue + '.allpairs.txt'
		f = open(file_name)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			gene_id = data[0]
			snp_id = data[1]
			if gene_id not in valid_genes or snp_id not in valid_variants:
				continue
			abs_dist = abs(float(data[2]))
			if abs_dist > tss_distance_thresh:
				continue
			test_id = gene_id + ':' + snp_id
			if test_id not in tests:
				tests[test_id] = 1.0
			else:
				tests[test_id] = tests[test_id] + 1.0
		f.close()
	t = open(all_tests_file, 'w')
	t.write('gene_id\tvariant_id\n')
	for test_id in tests.keys():
		if tests[test_id] != len(tissues):
			continue
		test_info = test_id.split(':')
		gene_id = test_info[0]
		variant_id = test_info[1]
		t.write(gene_id + '\t' + variant_id + '\n')
	t.close()


def print_sample_surveyed_covariates(sample_names, gtex_individual_information_file, output_file):
	f = open(gtex_individual_information_file)
	column_values = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if len(data) != 175:
			continue
		if head_count == 0:
			head_count = head_count + 1
			for i, ele in enumerate(data):
				column_values.append({})
			continue
		for i, ele in enumerate(data):
			if ele not in column_values[i]:
				column_values[i][ele] = 0
			column_values[i][ele] = column_values[i][ele] + 1
	f.close()
	valid_columns = {}
	for itera in range(len(column_values)):
		valid = True
		values = column_values[itera]
		if len(values) == 1 or len(values) > 5:
			valid = False
		for value in values.keys():
			if value != '1' and value != '0' and value != '96' and value != '97' and value != '98' and value != '99':
				valid = False
		if valid == True:
			# Limit to at least 5 positive examples
			if '1' in values and values['1'] > 4:
				valid_columns[itera] = 1
	dicti = {}
	head_count = 0
	f = open(gtex_individual_information_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			if len(data) != 175: 
				print('assumption error')
			header = []
			for i, ele in enumerate(data):
				if i in valid_columns:
					header.append(ele)
			header = np.asarray(header)
			continue
		indi_id = data[0]
		if len(data) != 175:
			vec = np.asarray(['99']*len(header))
			dicti[indi_id] = vec
			continue
		vec = []
		for i, ele in enumerate(data):
			if i in valid_columns:
				vec.append(ele)
		vec = np.asarray(vec)
		if len(vec) != len(header):
			print('assumptino eroror')
			pdb.set_trace()
		dicti[indi_id] = vec
	f.close()
	t = open(output_file,'w')
	t.write('sample_id\t' + '\t'.join(header) + '\n')
	for sample_id in sample_names:
		indi = sample_id.split(':')[0]
		t.write(sample_id + '\t' + '\t'.join(dicti[indi]) + '\n')
	t.close()

def create_sample_name_mapping(sample_names, tissues, tissues_alt, gtex_tpm_dir):
	mapping = {}
	reverse_mapping = {}
	for sample_name in sample_names:
		mapping[sample_name] = 'Null'
	# Loop through each tissue and fill in the tpm matrix
	for tissue_index, tissue in enumerate(tissues):
		tissue_alt = tissues_alt[tissue_index]
		# Stream tpm file in this tissue
		tpm_file = gtex_tpm_dir + tissue_alt + '.txt'
		head_count = 0
		f = open(tpm_file)
		for line in f:
			line = line.rstrip()
			data = line.split()
			# Header
			if head_count == 0:
				head_count = head_count + 1
				# Get ordered list of sample names
				sample_names = []
				for sample_name_temp in data[2:]:
					our_sample_name = sample_name_temp.split('-')[0] + '-' + sample_name_temp.split('-')[1] + ':' + tissue
					if our_sample_name in mapping:
						if mapping[our_sample_name] != 'Null':
							print('assumption erorro')
							pdb.set_trace()
						mapping[our_sample_name] = sample_name_temp
						reverse_mapping[sample_name_temp] = our_sample_name
					#sample_names.append(sample_name)
				continue
			break
		f.close()
	if len(mapping) != len(reverse_mapping):
		print('assumption erororr')
		pdb.set_trace()
	for sample_name in sample_names:
		if mapping[sample_name] == 'Null':
			print('assumption erororroororororo')
			pdb.set_trace()
	return reverse_mapping


def print_sample_technical_covariates(tissues, tissues_alt, sample_names, gtex_tpm_dir, gtex_sample_information_file, output_file):
	sample_name_mapping = create_sample_name_mapping(sample_names, tissues, tissues_alt, gtex_tpm_dir)
	sample_dicti = {}
	for sample in sample_names:
		sample_dicti[sample] = 0
	f = open(gtex_sample_information_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count =head_count + 1
			header = np.copy(data)
			continue
		sample_id = data[0]
		if len(data) != len(header):
			continue
		if sample_id not in sample_name_mapping:
			continue
		our_sample_id = sample_name_mapping[sample_id]
		sample_dicti[our_sample_id] = 1
	f.close()
	for key_name in sample_dicti.keys():
		if sample_dicti[key_name] != 1:
			print('assumption error')
			pdb.set_trace()

	covariate_names = ['SMNABTCH', 'SMRIN', 'SMGEBTCH', 'SMNTERRT', 'SMGNSDTC', 'SMSPLTRD', 'SMRRNART', 'SMEXNCRT', 'SMMPPD']
	num_cov = len(covariate_names)
	sample_to_cov = {}
	for sample in sample_names:
		sample_to_cov[sample] = np.zeros(num_cov)
	f = open(gtex_sample_information_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = np.asarray(line.split('\t'))
		if head_count == 0:
			head_count = head_count + 1
			header = np.copy(data)
			cov_indices = []
			for covariate_name in covariate_names:
				temp = np.where(header==covariate_name)[0]
				if len(temp) != 1:
					print('assumption eroror')
				cov_indices.append(temp[0])
			cov_indices = np.asarray(cov_indices)
			if len(cov_indices) != num_cov:
				print('assumption error')
				pdb.set_trace()
			continue
		sample_id = data[0]
		if sample_id not in sample_name_mapping:
			continue
		our_sample_id = sample_name_mapping[sample_id]
		sample_to_cov[our_sample_id] = data[cov_indices]
	f.close()
	t = open(output_file, 'w')
	t.write('sample_id\t' + '\t'.join(covariate_names) + '\n')
	for sample in sample_names:
		t.write(sample + '\t' + '\t'.join(sample_to_cov[sample]) + '\n')
	t.close()


def extract_categorical_feature(raw_feature_cov):
	all_features = []
	for raw_feature_str in raw_feature_cov:
		for raw_feature in raw_feature_str.split(','):
			all_features.append(raw_feature)
	unique_categories = np.unique(np.asarray(all_features))
	mapping = {}
	for i, category in enumerate(unique_categories):
		mapping[category] = i

	cat_feature = np.zeros((len(raw_feature_cov), len(unique_categories) -1))

	for i, raw_feature_string in enumerate(raw_feature_cov):
		for raw_feature in raw_feature_string.split(','):
			column_num = mapping[raw_feature]
			if column_num == (len(unique_categories) -1):
				continue
			cat_feature[i, column_num] = 1.0
	return(cat_feature)


def extract_technical_covariates(technical_covariate_file):
	raw_cov = np.loadtxt(technical_covariate_file, dtype=str, delimiter='\t')
	header = raw_cov[0,1:]
	raw_cov = raw_cov[1:, 1:]
	# a bit manual here
	feature_types = ['skip', 'real', 'cat', 'skip', 'skip', 'skip', 'skip', 'skip', 'real']
	# error checking
	if len(feature_types) != raw_cov.shape[1]:
		print('assumption error')
	feature_arr = []
	for feature_iter in range(len(feature_types)):
		feature_type = feature_types[feature_iter]
		raw_feature_cov = raw_cov[:, feature_iter]
		if feature_type == 'cat':
			cat_feature = extract_categorical_feature(raw_feature_cov)
			feature_arr.append(cat_feature)
		elif feature_type == 'real':
			feature_arr.append(raw_feature_cov[:, np.newaxis].astype(float))
		elif feature_type == 'skip':
			continue
		else:
			print('Error: shouldnt be here')
			pdb.set_trace()
	return np.hstack(feature_arr)



def extract_sample_covariates(covariate_file):
	raw_cov = np.loadtxt(covariate_file, dtype=str, delimiter='\t')
	header = raw_cov[0,1:]
	raw_cov = raw_cov[1:, 1:]
	# a bit manual here
	feature_types = ['skip', 'cat', 'skip', 'skip', 'skip', 'skip', 'skip', 'skip', 'real', 'skip', 'skip', 'skip', 'skip', 'skip', 'skip', 'skip']
	# error checking
	if len(feature_types) != raw_cov.shape[1]:
		print('assumption error')
	feature_arr = []
	for feature_iter in range(len(feature_types)):
		feature_type = feature_types[feature_iter]
		raw_feature_cov = raw_cov[:, feature_iter]
		if feature_type == 'cat':
			cat_feature = extract_categorical_feature(raw_feature_cov)
			feature_arr.append(cat_feature)
		elif feature_type == 'real':
			feature_arr.append(raw_feature_cov[:, np.newaxis].astype(float))
		elif feature_type == 'skip':
			continue
		else:
			print('Error: shouldnt be here')
			pdb.set_trace()
	return np.hstack(feature_arr)

def regress_out_technical_covariates_from_gene_expression(input_expression_file, output_residual_expression_file, sample_covariate_file, technical_covariate_file):
	# A bit of manual covariate extraction.. sue me
	technical_covariates = extract_technical_covariates(technical_covariate_file)
	donor_covariates = extract_sample_covariates(sample_covariate_file)
	covariates = np.hstack((donor_covariates, technical_covariates))

	expr_raw = np.loadtxt(input_expression_file, dtype=str, delimiter='\t')
	expr = expr_raw[1:,1:].astype(float)
	reg = LinearRegression().fit(covariates,np.transpose(expr))
	prediction = np.transpose(reg.predict(covariates))
	residual = expr - prediction

	f = open(input_expression_file)
	t = open(output_residual_expression_file,'w')
	counter = 0
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		ensamble_id = data[0]
		t.write(ensamble_id + '\t' + '\t'.join(residual[counter, :].astype(str)) + '\n')
		counter = counter + 1
	f.close()
	t.close()

def get_xcell_covariates(sample_names, tissues, tissues_alt, gtex_tpm_dir, sample_to_index, gtex_xcell_enrichment_file):
	sample_name_mapping = create_sample_name_mapping(sample_names, tissues, tissues_alt, gtex_tpm_dir)
	xcell_covariate_dicti = {}
	xcell_data = np.loadtxt(gtex_xcell_enrichment_file, dtype=str,delimiter='\t')
	xcell_cell_types = xcell_data[1:,0]
	xcell_sample_names = xcell_data[0,1:]
	xcell_enrichments = xcell_data[1:, 1:]
	num_samples = len(xcell_sample_names)
	for sample_num in range(num_samples):
		xcell_sample_name = xcell_sample_names[sample_num]
		if xcell_sample_name not in sample_name_mapping:
			continue
		our_sample_name = sample_name_mapping[xcell_sample_name]
		if our_sample_name in xcell_covariate_dicti:
			print('assumption erroror')
			pdb.set_trace()
		xcell_covariate_dicti[our_sample_name] = xcell_enrichments[:, sample_num]
	if len(sample_to_index) != len(xcell_covariate_dicti):
		print('assumption erororo!')
		pdb.set_trace()
	return xcell_covariate_dicti, xcell_cell_types

tissues_file = sys.argv[1]
gtex_expression_dir = sys.argv[2]
gtex_tpm_dir = sys.argv[3]
gtex_covariate_dir = sys.argv[4]
gtex_genotype_dir = sys.argv[5]
gtex_egene_dir = sys.argv[6]
gtex_individual_information_file = sys.argv[7]
gtex_sample_information_file = sys.argv[8]
gtex_eqtl_dir = sys.argv[9]
gtex_xcell_enrichment_file = sys.argv[10]
processed_data_dir = sys.argv[11]
tissues, tissues_alt = get_tissues(tissues_file)


# Extract file of sample names
sample_name_file = processed_data_dir + 'sample_names.txt'
sample_names, sample_to_index, individual_covariates = get_sample_names(tissues, gtex_expression_dir, sample_name_file, gtex_individual_information_file)

#print_individual_id_file_and_tissue_file(sample_names, processed_data_dir + 'individual_id.txt', processed_data_dir + 'sample_tissue_names.txt')

xcell_covariates, xcell_covariate_names = get_xcell_covariates(sample_names, tissues, tissues_alt, gtex_tpm_dir, sample_to_index, gtex_xcell_enrichment_file)

sample_covariate_file = processed_data_dir + 'sample_covariates.txt'
print_sample_covariates(sample_names, individual_covariates, xcell_covariates, xcell_covariate_names, sample_covariate_file)

#print_sample_surveyed_covariates(sample_names, gtex_individual_information_file, processed_data_dir + 'sample_surveyed_covariates.txt')

technical_covariate_file = processed_data_dir + 'sample_technical_covariates.txt'
#print_sample_technical_covariates(tissues, tissues_alt, sample_names, gtex_tpm_dir, gtex_sample_information_file, technical_covariate_file)

# We are going to limit analysis to genes tested in all tissues
#genes_tested_in_all_tissues = get_genes_tested_in_all_tissues(tissues, gtex_expression_dir)


# Generate TPM expression matrix
tpm_expression_matrix_file = processed_data_dir + 'cross_tissue_tpm.txt'
#generate_tpm_expression_matrix(tissues, tissues_alt, sample_names, genes_tested_in_all_tissues, gtex_tpm_dir, tpm_expression_matrix_file)

# Quantile normalize and standardize TPM expression matrix
standardized_tpm_expression_matrix_file = processed_data_dir + 'cross_tissue_tpm_standardized.txt'
#standardize_expression(tpm_expression_matrix_file, standardized_tpm_expression_matrix_file)


# Extract covariates (expression pcs)
num_expression_pcs = 50
covariate_file = processed_data_dir + 'covariates.txt'
#extract_covariates(standardized_tpm_expression_matrix_file, covariate_file, num_expression_pcs, gtex_covariate_dir, sample_names, tissues)



# Limit to genes in our analysis
#valid_genes = get_genes_we_have_expression_data_for(standardized_tpm_expression_matrix_file)

# Limit to variants we have genoytpe data for
#valid_variants = get_variants_we_have_genotype_data_for(gtex_genotype_dir)

tss_distance_thresh=50000.0
all_tests_file = processed_data_dir + 'all_tests_unordered.txt'
#extract_all_variant_gene_pairs(tissues, gtex_eqtl_dir, tss_distance_thresh, valid_genes, valid_variants, all_tests_file)




# Regress out technical covariates from gene expression
standardized_tpm_expression_technical_cov_residual_matrix_file = processed_data_dir + 'cross_tissue_tpm_technical_covariate_residuals.txt'
#regress_out_technical_covariates_from_gene_expression(standardized_tpm_expression_matrix_file, standardized_tpm_expression_technical_cov_residual_matrix_file, sample_covariate_file, technical_covariate_file)



# Extract covariates (expression pcs)
num_expression_pcs = 80
covariate_file = processed_data_dir + 'residual_expression_covariates.txt'
#extract_covariates(standardized_tpm_expression_technical_cov_residual_matrix_file, covariate_file, num_expression_pcs, gtex_covariate_dir, sample_names, tissues)














