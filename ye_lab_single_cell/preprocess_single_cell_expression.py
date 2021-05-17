import numpy as np 
import os
import sys
import pdb
import h5py
from sklearn.decomposition import PCA
import scanpy as sc
from anndata import AnnData
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.linear_model import LinearRegression
from sklearn.cluster import KMeans
import rnaseqnorm
import pandas as pd



def extract_protein_coding_known_autosomal_genes(gene_struct, gene_annotation_file):
	'''
	# Simple error checking to understand gene struct
	col_names = gene_struct.columns
	ensamble_ids = gene_struct[col_names[0]]
	for col_name in col_names:
		if len(np.unique(gene_struct[col_name])) != 32738:
			print('assumption errror')
			pdb.set_trace()
		if np.array_equal(gene_struct[col_name], ensamble_ids) == False:
			print('assumption errorr')
			pdb.set_trace()
	'''
	# Make dictionary list of genes that are protein-coding, known, and autosomal
	valid_genes = {}
	f = open(gene_annotation_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Simple error check
		if len(data) != 8 and len(data) != 9:
			print('assumption errror in processing gene annotation file')
			pdb.set_trace()
		gene_id = data[5]
		# Some more error checks
		start = int(data[2])
		end = int(data[3])
		if start > end:
			print('assumption errror in processing gene annotation file')
			pdb.set_trace()
		if data[6] != 'protein_coding' or data[7] != 'KNOWN':
			continue
		if data[1] == 'chrX' or data[1] == 'chrY' or data[1] == 'chrM':
			continue
		valid_genes[data[5] + '_' + data[0].split('.')[0]] = 1
	f.close()
	# Make binary vectory corresponding to ordered list of whehter our genes are protein-coding, autosomal and known
	binary_vector = []
	ensamble_ids = gene_struct[gene_struct.columns[0]]
	gene_ids = gene_struct.index
	num_genes = len(gene_ids)
	for gene_num in range(num_genes):
		ensamble_id = ensamble_ids[gene_num]
		gene_id = gene_ids[gene_num]
		if gene_id + '_' + ensamble_id not in valid_genes:
			binary_vector.append(False)
		else:
			binary_vector.append(True)
	return np.asarray(binary_vector)




# Generate expression PC loadings and variance explained of those expression PCs
def generate_pca_scores_and_variance_explained(filtered_standardized_sc_expression_file, num_pcs, filtered_cells_pca_file, filtered_cells_pca_ve_file):
	# Load in data
	X = np.loadtxt(filtered_standardized_sc_expression_file)

	# Run PCA (via SVD)
	#uuu, sss, vh = np.linalg.svd(np.transpose(X), full_matrices=False)
	#svd_loadings = np.transpose(vh)[:,:num_pcs]
	#ve = (np.square(sss)/np.sum(np.square(sss)))[:num_pcs]

	# Faster in sklearn
	_pca = PCA(n_components=num_pcs, svd_solver='arpack')
	svd_loadings = _pca.fit_transform(X)
	ve = _pca.explained_variance_ratio_

	# Save to output file
	np.savetxt(filtered_cells_pca_file, svd_loadings, fmt="%s", delimiter='\t')

	# Compute variance explained
	np.savetxt(filtered_cells_pca_ve_file, ve, fmt="%s", delimiter='\n')

# Generate pseudobulk covariate file (along with ordered list of pseudobulk samples)
def generate_pseudobulk_covariate_file(adata, pseudobulk_covariate_file):
	# Loop through cell covariate file to create  list of individual-cell_type pairs (these are our samples)
	samples = {}
	num_cells = adata.obs.shape[0]
	sample_counts = {}
	for cell_num in range(num_cells):
		# Extract relevent fields
		cell_type = adata.obs.ct_cov[cell_num]
		race = adata.obs.pop_cov[cell_num]
		ind_id = adata.obs.ind_cov[cell_num]
		# Come up with sample id name
		cell_type = '_'.join(cell_type.split(' '))
		sample_name = ind_id + ':' + cell_type
		if sample_name not in samples:  # Add sample to list
			samples[sample_name] = race
		else:  # Simple error checking making sure ever individual was mapped to the same race
			if samples[sample_name] != race:
				print('assumption errro')
		if sample_name not in sample_counts:
			sample_counts[sample_name] = 0
		sample_counts[sample_name] = sample_counts[sample_name] + 1
	# Print to output file
	t = open(pseudobulk_covariate_file, 'w')
	t.write('sample_name\tct_cov\tpop_cov\tind_cov\tnum_cells\tct_cov_readable\n')
	# Loop through samples and print to output file
	for sample_name in samples.keys():
		# Extract relevent fields
		ind_cov = sample_name.split(':')[0]
		ct_cov = ' '.join(sample_name.split(':')[1].split('_'))
		pop_cov = samples[sample_name]
		# Print to output file
		t.write(sample_name + '\t' + ct_cov + '\t' + pop_cov + '\t' + ind_cov + '\t' + str(sample_counts[sample_name]) + '\t' + sample_name.split(':')[1] + '\n')
	t.close()

def generate_pseudobulk_per_individual_covariate_file(adata, pseudobulk_covariate_file):
	# Loop through cell covariate file to create  list of individuals (these are our samples)
	samples = {}
	num_cells = adata.obs.shape[0]
	sample_counts = {}
	for cell_num in range(num_cells):
		# Extract relevent fields
		cell_type = adata.obs.ct_cov[cell_num]
		race = adata.obs.pop_cov[cell_num]
		ind_id = adata.obs.ind_cov[cell_num]
		# Come up with sample id name
		cell_type = '_'.join(cell_type.split(' '))
		sample_name = ind_id
		if sample_name not in samples:  # Add sample to list
			samples[sample_name] = race
		else:  # Simple error checking making sure ever individual was mapped to the same race
			if samples[sample_name] != race:
				print('assumption errro')
		if sample_name not in sample_counts:
			sample_counts[sample_name] = 0
		sample_counts[sample_name] = sample_counts[sample_name] + 1
	# Print to output file
	t = open(pseudobulk_covariate_file, 'w')
	t.write('sample_name\tpop_cov\tind_cov\tnum_cells\n')
	# Loop through samples and print to output file
	for sample_name in samples.keys():
		# Extract relevent fields
		#ind_cov = sample_name.split(':')[0]
		#ct_cov = ' '.join(sample_name.split(':')[1].split('_'))
		pop_cov = samples[sample_name]
		# Print to output file
		t.write(sample_name + '\t' + pop_cov + '\t' + sample_name + '\t' + str(sample_counts[sample_name]) + '\n')
	t.close()

# Generate summed pseudobulk counts
def generate_pseudobulk_counts(adata, raw_pseudobulk_expression_file, pseudobulk_covariate_file):
	# Goal: convert adata.raw.X (cellsXgenes) to pseudobulk_matrix (pseudobulk_samplesXgenes)
	# First create mapping from pseudobulk sample name to index. Also count number of pseudobulk samples
	pseudobulk_sample_to_index = {}
	f = open(pseudobulk_covariate_file)
	head_count = 0
	num_pseudobulk_samples = 0

	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		sample_id = data[0]
		pseudobulk_sample_to_index[sample_id] = num_pseudobulk_samples
		num_pseudobulk_samples = num_pseudobulk_samples + 1
	f.close()
	# Now create array of length number of cells where each element is mapping from cell to pseudobulk sample
	cell_to_pseudobulk_sample_index = []
	num_cells = adata.obs.shape[0]
	for cell_num in range(num_cells):
		# Extract relevent fields
		cell_type = adata.obs.ct_cov[cell_num]
		ind_id = adata.obs.ind_cov[cell_num]
		# Come up with sample id name for this cll
		cell_type = '_'.join(cell_type.split(' '))
		sample_name = ind_id + ':' + cell_type
		# Map cells sample id to an index
		pseudobulk_index = pseudobulk_sample_to_index[sample_name]
		cell_to_pseudobulk_sample_index.append(pseudobulk_index)
	cell_to_pseudobulk_sample_index = np.asarray(cell_to_pseudobulk_sample_index)
	# Now generate and fill in pseudobulk counts matrix
	num_genes = adata.raw.X.shape[1]
	pseudobulk_matrix = np.zeros((num_pseudobulk_samples, num_genes))
	# convert raw cell counts into np matrix
	raw_cell_counts = adata.raw.X.toarray()
	# Loop through cells and genes and add cell-gene count to pseudobulk count matrix
	for cell_num in range(num_cells):
		pseudobulk_sample_index = cell_to_pseudobulk_sample_index[cell_num]
		pseudobulk_matrix[pseudobulk_sample_index, :] = pseudobulk_matrix[pseudobulk_sample_index, :] + raw_cell_counts[cell_num, :]
	np.savetxt(raw_pseudobulk_expression_file, pseudobulk_matrix, fmt="%s", delimiter='\t')

def generate_pseudobulk_per_individual_counts(adata, raw_pseudobulk_expression_file, pseudobulk_covariate_file):
	# Goal: convert adata.raw.X (cellsXgenes) to pseudobulk_matrix (pseudobulk_samplesXgenes)
	# First create mapping from pseudobulk sample name to index. Also count number of pseudobulk samples
	pseudobulk_sample_to_index = {}
	f = open(pseudobulk_covariate_file)
	head_count = 0
	num_pseudobulk_samples = 0

	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		sample_id = data[0]
		pseudobulk_sample_to_index[sample_id] = num_pseudobulk_samples
		num_pseudobulk_samples = num_pseudobulk_samples + 1
	f.close()
	# Now create array of length number of cells where each element is mapping from cell to pseudobulk sample
	cell_to_pseudobulk_sample_index = []
	num_cells = adata.obs.shape[0]
	for cell_num in range(num_cells):
		# Extract relevent fields
		#cell_type = adata.obs.ct_cov[cell_num]
		ind_id = adata.obs.ind_cov[cell_num]
		# Come up with sample id name for this cll
		#cell_type = '_'.join(cell_type.split(' '))
		sample_name = ind_id
		# Map cells sample id to an index
		pseudobulk_index = pseudobulk_sample_to_index[sample_name]
		cell_to_pseudobulk_sample_index.append(pseudobulk_index)
	cell_to_pseudobulk_sample_index = np.asarray(cell_to_pseudobulk_sample_index)
	# Now generate and fill in pseudobulk counts matrix
	num_genes = adata.raw.X.shape[1]
	pseudobulk_matrix = np.zeros((num_pseudobulk_samples, num_genes))
	# convert raw cell counts into np matrix
	raw_cell_counts = adata.raw.X.toarray()
	# Loop through cells and genes and add cell-gene count to pseudobulk count matrix
	for cell_num in range(num_cells):
		pseudobulk_sample_index = cell_to_pseudobulk_sample_index[cell_num]
		pseudobulk_matrix[pseudobulk_sample_index, :] = pseudobulk_matrix[pseudobulk_sample_index, :] + raw_cell_counts[cell_num, :]
	np.savetxt(raw_pseudobulk_expression_file, pseudobulk_matrix, fmt="%s", delimiter='\t')


# Standardize summed pseudobulk counts (samples X genes)
def standardize_pseudobulk_counts(raw_pseudobulk_expression_file, standardized_pseudobulk_expression_file):
	# Load in raw pseodobulk expression data
	raw_pseudobulk_expression = np.loadtxt(raw_pseudobulk_expression_file)
	num_samples, num_genes = raw_pseudobulk_expression.shape

	# initialize and fill in normalized pseudobulk expression data
	normalized_pseudobulk_expression = np.zeros((num_samples, num_genes))
	for sample_num in range(num_samples):
		normalized_pseudobulk_expression[sample_num, :] = np.log(10000.0*(raw_pseudobulk_expression[sample_num,:]/np.sum(raw_pseudobulk_expression[sample_num,:])) + 1.0)
	
	# initialize and fill in standardized pseudobulk expression data
	standardized_pseudobulk_expression = np.zeros((num_samples, num_genes))
	for gene_num in range(num_genes):
		standardized_pseudobulk_expression[:, gene_num] = (normalized_pseudobulk_expression[:, gene_num] - np.mean(normalized_pseudobulk_expression[:, gene_num]))/np.std(normalized_pseudobulk_expression[:, gene_num])

	np.savetxt(standardized_pseudobulk_expression_file, np.round(standardized_pseudobulk_expression, decimals=5), fmt="%s", delimiter='\t')


def generate_pseudobulk_expression_data_wrapper(adata, raw_pseudobulk_expression_file, standardized_pseudobulk_expression_file, pseudobulk_covariate_file):
	# Generate pseudobulk covariate file (along with ordered list of pseudobulk samples)
	generate_pseudobulk_covariate_file(adata, pseudobulk_covariate_file)

	# Generate summed pseudobulk counts
	generate_pseudobulk_counts(adata, raw_pseudobulk_expression_file, pseudobulk_covariate_file)

	# Standardize summed pseudobulk counts (samples X genes)
	standardize_pseudobulk_counts(raw_pseudobulk_expression_file, standardized_pseudobulk_expression_file)

def generate_pseudobulk_per_individual_expression_data_wrapper(adata, raw_pseudobulk_expression_file, standardized_pseudobulk_expression_file, pseudobulk_covariate_file):
	# Generate pseudobulk covariate file (along with ordered list of pseudobulk samples)
	generate_pseudobulk_per_individual_covariate_file(adata, pseudobulk_covariate_file)

	# Generate summed pseudobulk counts
	generate_pseudobulk_per_individual_counts(adata, raw_pseudobulk_expression_file, pseudobulk_covariate_file)

	# Standardize summed pseudobulk counts (samples X genes)
	standardize_pseudobulk_counts(raw_pseudobulk_expression_file, standardized_pseudobulk_expression_file)

def center_each_column(raw_counts):
	num_genes = raw_counts.shape[1]
	num_cells = raw_counts.shape[0]
	new_mat = np.zeros((num_cells, num_genes))
	for gene_num in range(num_genes):
		new_mat[:, gene_num] = (raw_counts[:, gene_num] - np.mean(raw_counts[:, gene_num]))
	return new_mat


def generate_cell_type_sc_expression_data_wrapper(cell_type, adata, standardized_cell_type_expression_file, standardized_cell_type_covariate_file, raw_cell_type_expression_file, centered_raw_cell_type_expression_file):
	############################
	# Filter to cells of this cell type
	###############################
	ct_data = adata[adata.obs.ct_cov == cell_type, :]

	############################
	# Save covariates to output file
	###############################
	np.savetxt(standardized_cell_type_covariate_file, ct_data.obs, fmt="%s", delimiter='\t', header='\t'.join(ct_data.obs.columns), comments='')

	############################
	# Save raw expression to output file
	###############################
	np.savetxt(raw_cell_type_expression_file, ct_data.raw.X.toarray(), fmt="%s", delimiter='\t')

	centered_raw_counts = center_each_column(ct_data.raw.X.toarray())
	np.savetxt(centered_raw_cell_type_expression_file, centered_raw_counts, fmt="%s", delimiter='\t')
	
	######################	
	# Normalize data
	#######################
	ct_data_clean = AnnData(ct_data.raw.X.toarray())
	sc.pp.normalize_total(ct_data_clean, target_sum=1e4)
	sc.pp.log1p(ct_data_clean)
	sc.pp.scale(ct_data_clean, max_value=1000)
	############################
	# Save expression to output
	###########################
	np.savetxt(standardized_cell_type_expression_file, ct_data_clean.X, fmt="%s", delimiter='\t')

def calculate_deviance_residuals(adata):
	adata.X = adata.X.toarray()
	num_cells, num_genes = adata.X.shape
	total_cell_counts = np.sum(adata.X,axis=1)
	total_gene_counts = np.sum(adata.X, axis=0)
	gene_frac = total_gene_counts/np.sum(total_gene_counts)


	'''
	for j in range(num_genes):
		# Observed count for gene, individual pair
		y_j = adata.X[:,j]
		# Compute predicted fraction assigned to gene (according to MLE)
		pi_j = gene_frac[j]
		# Compute predicted fraction assigned to gene, individual pair
		mu_j = total_cell_counts*pi_j
		# Compute deviance residual
		term_1 = 2.0*y_j*np.log(y_j/mu_j) 
		term_1[np.isnan(term_1)] = 0.0
		term_2 = 2.0*(total_cell_counts - y_j)*np.log((total_cell_counts - y_j)/(total_cell_counts - mu_j))
		deviance_residual_j = np.sign(y_j - mu_j)*np.sqrt(term_1 + term_2)
	'''
	# Doing pearson resid instead
	for j in range(num_genes):
		# Observed count for gene, individual pair
		y_j = adata.X[:,j]
		# Compute predicted fraction assigned to gene (according to MLE)
		pi_j = gene_frac[j]

		# Compute predicted fraction assigned to gene, individual pair
		mu_j = total_cell_counts*pi_j

		# Compute pearson resid
		numerator_term = y_j - mu_j

		denomenator_term= mu_j - ((1.0/total_cell_counts)*np.square(mu_j))

		pearson_resid_j = numerator_term/np.sqrt(denomenator_term)
		adata.X[:, j] = pearson_resid_j
		standardized_pearson_resid = (pearson_resid_j - np.mean(pearson_resid_j))/np.std(pearson_resid_j)
		adata.X[:, j] = standardized_pearson_resid
	return adata

def create_cell_similarity_matrix_with_pearson_correlation(individual_expression):
	# Compute pearson correlation of each cell 
	pearson_corr_mat = np.corrcoef(individual_expression)
	return pearson_corr_mat

def create_cell_similarity_matrix_with_euclidean_distance(individual_expression):
	num_samples = individual_expression.shape[0]
	distance_mat = euclidean_distances(individual_expression, individual_expression)
	'''
	for i in range(num_samples):
		print(i)
		for j in range(num_samples):
			dist = np.linalg.norm(individual_expression[i, :] - individual_expression[j, :])
			distance_mat[i, j] = dist
	'''
	return -distance_mat

def standardize_columns(expr_mat):
	standardized_expr_mat = np.zeros(expr_mat.shape)
	num_cols = expr_mat.shape[1]
	for col_num in range(num_cols):
		standardized_expr_mat[:, col_num] = (expr_mat[:, col_num] - np.mean(expr_mat[:, col_num]))/np.std(expr_mat[:, col_num])
	return standardized_expr_mat

def create_knn_mapping(adata, k, knn_method):
	# Initialize mapping from cell id to column positions that are the kNN to that cell id
	knn_mapping_to_indices = {}
	# Initialize mapping from cell id to cell_ids that are the kNN to that cell id
	knn_mapping_to_cell_ids = {}
	# Knn mapping done seperately for each individual
	cell_individuals = np.asarray(adata.obs['ind_cov'])
	unique_individuals = np.unique(cell_individuals)
	cell_names = np.asarray(adata.obs['cell_id'])
	cell_names_to_index = {}
	for i, cell_name in enumerate(cell_names):
		cell_names_to_index[cell_name] = i

	standardized_pca_loadings = standardize_columns(adata.obsm['X_pca'])

	# Loop through individuals
	for individual in unique_individuals:
		print(individual)
		# Get indices of cells corresponding to this individual
		individual_indices = np.where(cell_individuals == individual)[0]
		# Simple error checking
		if 2.0*k > len(individual_indices):
			print('assumption error: not enough cells per individual')
			pdb.set_trace()
		# Filter data to strictly cells from this individual
		individual_expression = adata.X[individual_indices, :]
		individual_pca_loadings = standardized_pca_loadings[individual_indices, :20]
		individual_cell_names = cell_names[individual_indices]
		# Simple error checking
		if len(np.unique(cell_individuals[individual_indices])) != 1:
			print('assumption error')
			pdb.set_trace()
		# compute similarity between cells
		if knn_method == 'pearson_full':
			cell_similarity_matrix = create_cell_similarity_matrix_with_pearson_correlation(individual_expression)
		elif knn_method == 'euclidean_pca':
			cell_similarity_matrix = create_cell_similarity_matrix_with_euclidean_distance(individual_pca_loadings)
		# Loop through cells
		for i, cell_name in enumerate(individual_cell_names):
			cells_similarity_vector = cell_similarity_matrix[i,:]
			top_k_cell_indices = cells_similarity_vector.argsort()[-(k+1):][::-1]
			if len(top_k_cell_indices) != (k+1):
				print('assumption error')
				pdb.set_trace()
			top_k_cell_names = []
			top_k_indices = []
			for cell_index in top_k_cell_indices:
				if cell_index == i:
					continue
				top_k_cell_names.append(individual_cell_names[cell_index])
				index = cell_names_to_index[individual_cell_names[cell_index]]
				#index = np.where(cell_names == individual_cell_names[cell_index])[0]
				#if len(index) != 1:
					#print('assumption error')
					#pdb.set_trace()
				top_k_indices.append(index)
			top_k_cell_names = np.asarray(top_k_cell_names)
			top_k_indices = np.asarray(top_k_indices)
			# Error checking
			if len(top_k_indices) != k or len(np.unique(top_k_indices)) != k:
				print('assumption error')
				pdb.set_trace()
			if len(top_k_cell_names) != k or len(np.unique(top_k_cell_names)) != k:
				print('assumption error')
				pdb.set_trace()
			# Add info to dictionary
			if cell_name in knn_mapping_to_indices:
				print('assumption error')
				pdb.set_trace()
			if cell_name in knn_mapping_to_cell_ids:
				print('assumption erorr')
				pdb.set_trace()
			knn_mapping_to_indices[cell_name] = top_k_indices
			knn_mapping_to_cell_ids[cell_name] = top_k_cell_names
	if len(cell_names) != len(knn_mapping_to_indices):
		print('aasssumptionerror')
		pdb.set_trace()
	return knn_mapping_to_indices, knn_mapping_to_cell_ids, cell_names

def print_knn_mapping(knn_mapping_to_indices, knn_mapping_to_cell_ids, ordered_cell_ids, knn_mapping_file):
	t = open(knn_mapping_file,'w')
	t.write('cell_id\tneighbor_indices\tneighbor_cell_ids\n')
	for cell_id in ordered_cell_ids:
		neighbor_indices = knn_mapping_to_indices[cell_id].astype(str)
		neighbor_cell_ids = knn_mapping_to_cell_ids[cell_id]
		t.write(cell_id + '\t' + ','.join(neighbor_indices) + '\t' + ','.join(neighbor_cell_ids) + '\n')
	t.close()

def print_distance_weighted_knn_mapping(knn_mapping_to_indices, knn_mapping_to_cell_ids, knn_mapping_to_weights, ordered_cell_ids, cell_id_to_kernel_width, knn_mapping_file):
	t = open(knn_mapping_file,'w')
	t.write('cell_id\tcell_id_kernel_width\tneighbor_indices\tneighbor_cell_ids\tneighbor_weights\n')
	for cell_id in ordered_cell_ids:
		kernel_width = cell_id_to_kernel_width[cell_id]
		neighbor_indices = knn_mapping_to_indices[cell_id].astype(str)
		neighbor_cell_ids = knn_mapping_to_cell_ids[cell_id]
		neighbor_weights = knn_mapping_to_weights[cell_id].astype(str)
		t.write(cell_id + '\t' + str(kernel_width) + '\t' + ','.join(neighbor_indices) + '\t' + ','.join(neighbor_cell_ids) + '\t' + ','.join(neighbor_weights) + '\n')
	t.close()


def temp_load_knn_mapping(knn_mapping_file):
	knn_mapping_to_indices = {}
	knn_mapping_to_cell_ids = {}
	ordered_cell_ids = []
	f = open(knn_mapping_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		cell_id = data[0]
		neighbor_indices = np.asarray(data[1].split(',')).astype(int)
		neighbor_cell_ids = np.asarray(data[2].split(','))
		ordered_cell_ids.append(cell_id)
		knn_mapping_to_indices[cell_id] = neighbor_indices
		knn_mapping_to_cell_ids[cell_id] = neighbor_cell_ids
	f.close()
	ordered_cell_ids = np.asarray(ordered_cell_ids)
	return knn_mapping_to_indices, knn_mapping_to_cell_ids, ordered_cell_ids

def temp_load_distance_weighted_knn_mapping(knn_mapping_file):
	knn_mapping_to_indices = {}
	knn_mapping_to_cell_ids = {}
	knn_mapping_to_weights = {}
	cell_id_to_kernel_width = {}
	ordered_cell_ids = []
	f = open(knn_mapping_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		cell_id = data[0]
		kernel_width = float(data[1])
		neighbor_indices = np.asarray(data[2].split(',')).astype(int)
		neighbor_cell_ids = np.asarray(data[3].split(','))
		neighbor_weights = np.asarray(data[4].split(',')).astype(float)
		ordered_cell_ids.append(cell_id)
		knn_mapping_to_indices[cell_id] = neighbor_indices
		knn_mapping_to_cell_ids[cell_id] = neighbor_cell_ids
		knn_mapping_to_weights[cell_id] = neighbor_weights
		cell_id_to_kernel_width[cell_id] = kernel_width
	f.close()
	ordered_cell_ids = np.asarray(ordered_cell_ids)
	return knn_mapping_to_indices, knn_mapping_to_cell_ids, knn_mapping_to_weights, ordered_cell_ids, cell_id_to_kernel_width

def generate_knn_boosted_counts(adata, knn_mapping_to_indices, ordered_cell_ids, raw_knn_boosted_expression_file):
	# convert raw cell counts into np matrix
	raw_cell_counts = adata.raw.X.toarray()
	# Initialize knn-boosted raw cell counts
	knn_boosted_raw_cell_counts = np.zeros(raw_cell_counts.shape)

	# Loop through each cells
	num_cells = len(ordered_cell_ids)
	for cell_num in range(num_cells):
		# Get indices corresponding to KNN of this cell
		cell_name = ordered_cell_ids[cell_num]
		knn_indices = knn_mapping_to_indices[cell_name]
		# Sum across cell and its knn
		knn_boosted_raw_cell_counts[cell_num, :] = np.sum(raw_cell_counts[knn_indices,:],axis=0) + raw_cell_counts[cell_num,:]
	np.savetxt(raw_knn_boosted_expression_file, knn_boosted_raw_cell_counts, fmt="%s", delimiter='\t')

def generate_distance_weighted_knn_boosted_counts(adata, knn_mapping_to_indices, ordered_cell_ids, knn_mapping_to_weights, raw_knn_boosted_expression_file):
	# convert raw cell counts into np matrix
	raw_cell_counts = adata.raw.X.toarray()
	# Initialize knn-boosted raw cell counts
	knn_boosted_raw_cell_counts = np.zeros(raw_cell_counts.shape)

	# Loop through each cells
	num_cells = len(ordered_cell_ids)
	for cell_num in range(num_cells):
		# Get indices corresponding to KNN of this cell
		cell_name = ordered_cell_ids[cell_num]
		knn_indices = knn_mapping_to_indices[cell_name]
		knn_weights = knn_mapping_to_weights[cell_name]
		weight_denom = np.sum(knn_weights) + 1.0 # 1.0 comes from cell itself (which has weight of 1)
		knn_boosted_raw_cell_counts[cell_num, :] = knn_boosted_raw_cell_counts[cell_num, :] + (1.0/weight_denom)*raw_cell_counts[cell_num,:]
		for i, knn_index in enumerate(knn_indices):
			weight = knn_weights[i]
			knn_boosted_raw_cell_counts[cell_num, :] = knn_boosted_raw_cell_counts[cell_num, :] + (weight/weight_denom)*raw_cell_counts[knn_index,:]
		# Sum across cell and its knn
		#knn_boosted_raw_cell_counts[cell_num, :] = np.sum(raw_cell_counts[knn_indices,:],axis=0) + raw_cell_counts[cell_num,:]
	np.savetxt(raw_knn_boosted_expression_file, knn_boosted_raw_cell_counts, fmt="%s", delimiter='\t')

def print_distance_weighted_knn_mapping_cell_type_summary(adata, knn_mapping_file, knn_mapping_ct_summary_file):
	ordered_cell_types = np.asarray(adata.obs['ct_cov'])
	ordered_cell_ids = np.asarray(adata.obs['cell_id'])
	# create mapping from cell id to cell type
	cell_id_to_cell_type = {}
	for i, cell_id in enumerate(ordered_cell_ids):
		cell_id_to_cell_type[cell_id] = ordered_cell_types[i]
	# Get unique list of cell types
	unique_cell_types = sorted(np.unique(ordered_cell_types))
	num_cell_types = len(unique_cell_types)
	cell_type_to_pos = {}
	for i, cell_type in enumerate(unique_cell_types):
		cell_type_to_pos[cell_type] = i
	# Create dictionary to keep track of cell type knn counts
	mapping = {}
	for cell_type in unique_cell_types:
		mapping[cell_type] = np.zeros(num_cell_types)

	# Loop through knn file 
	f = open(knn_mapping_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		cell_id = data[0]
		cell_type_main = cell_id_to_cell_type[cell_id]
		knn_cell_ids = data[3].split(',')
		knn_weights = np.asarray(data[4].split(',')).astype(float)
		knn_normalized_weights = knn_weights/np.sum(knn_weights)
		for i, knn_cell_id in enumerate(knn_cell_ids):
			knn_cell_type = cell_id_to_cell_type[knn_cell_id]
			knn_cell_pos = cell_type_to_pos[knn_cell_type]
			mapping[cell_type_main][knn_cell_pos] = mapping[cell_type_main][knn_cell_pos] + knn_normalized_weights[i]
	f.close()
	t = open(knn_mapping_ct_summary_file, 'w')
	t.write('cell_type\t' + '\t'.join(unique_cell_types) + '\n')
	for cell_type in unique_cell_types:
		t.write(cell_type + '\t')
		frac_ct = mapping[cell_type]/np.sum(mapping[cell_type])
		t.write('\t'.join(frac_ct.astype(str)) + '\n')
	t.close()

def print_knn_mapping_batch_summary(adata, knn_mapping_file, knn_mapping_ct_summary_file):
	ordered_cell_types = np.asarray(adata.obs['batch_cov'])
	ordered_cell_ids = np.asarray(adata.obs['cell_id'])
	# create mapping from cell id to cell type
	cell_id_to_cell_type = {}
	for i, cell_id in enumerate(ordered_cell_ids):
		cell_id_to_cell_type[cell_id] = ordered_cell_types[i]
	# Get unique list of cell types
	unique_cell_types = sorted(np.unique(ordered_cell_types))
	num_cell_types = len(unique_cell_types)
	cell_type_to_pos = {}
	for i, cell_type in enumerate(unique_cell_types):
		cell_type_to_pos[cell_type] = i
	# Create dictionary to keep track of cell type knn counts
	mapping = {}
	for cell_type in unique_cell_types:
		mapping[cell_type] = np.zeros(num_cell_types)

	# Loop through knn file 
	f = open(knn_mapping_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		cell_id = data[0]
		cell_type_main = cell_id_to_cell_type[cell_id]
		knn_cell_ids = data[2].split(',')
		for knn_cell_id in knn_cell_ids:
			knn_cell_type = cell_id_to_cell_type[knn_cell_id]
			knn_cell_pos = cell_type_to_pos[knn_cell_type]
			mapping[cell_type_main][knn_cell_pos] = mapping[cell_type_main][knn_cell_pos] + 1
	f.close()
	t = open(knn_mapping_ct_summary_file, 'w')
	t.write('batch\t' + '\t'.join(unique_cell_types) + '\n')
	for cell_type in unique_cell_types:
		t.write(cell_type + '\t')
		frac_ct = mapping[cell_type]/np.sum(mapping[cell_type])
		t.write('\t'.join(frac_ct.astype(str)) + '\n')
	t.close()

def print_knn_mapping_cell_type_summary(adata, knn_mapping_file, knn_mapping_ct_summary_file):
	ordered_cell_types = np.asarray(adata.obs['ct_cov'])
	ordered_cell_ids = np.asarray(adata.obs['cell_id'])
	# create mapping from cell id to cell type
	cell_id_to_cell_type = {}
	for i, cell_id in enumerate(ordered_cell_ids):
		cell_id_to_cell_type[cell_id] = ordered_cell_types[i]
	# Get unique list of cell types
	unique_cell_types = sorted(np.unique(ordered_cell_types))
	num_cell_types = len(unique_cell_types)
	cell_type_to_pos = {}
	for i, cell_type in enumerate(unique_cell_types):
		cell_type_to_pos[cell_type] = i
	# Create dictionary to keep track of cell type knn counts
	mapping = {}
	for cell_type in unique_cell_types:
		mapping[cell_type] = np.zeros(num_cell_types)

	# Loop through knn file 
	f = open(knn_mapping_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		cell_id = data[0]
		cell_type_main = cell_id_to_cell_type[cell_id]
		knn_cell_ids = data[2].split(',')
		for knn_cell_id in knn_cell_ids:
			knn_cell_type = cell_id_to_cell_type[knn_cell_id]
			knn_cell_pos = cell_type_to_pos[knn_cell_type]
			mapping[cell_type_main][knn_cell_pos] = mapping[cell_type_main][knn_cell_pos] + 1
	f.close()
	t = open(knn_mapping_ct_summary_file, 'w')
	t.write('cell_type\t' + '\t'.join(unique_cell_types) + '\n')
	for cell_type in unique_cell_types:
		t.write(cell_type + '\t')
		frac_ct = mapping[cell_type]/np.sum(mapping[cell_type])
		t.write('\t'.join(frac_ct.astype(str)) + '\n')
	t.close()

def create_distance_weighted_knn_mapping(adata, k, knn_method):
	#Initialize mapping from cell id to kernell width
	cell_id_to_kernel_width = {}
	# Initialize mapping from cell id to column positions that are the kNN to that cell id
	knn_mapping_to_indices = {}
	# Initialize mapping from cell id to cell_ids that are the kNN to that cell id
	knn_mapping_to_cell_ids = {}
	# Initialize mapping from cell id weights of cell ids
	knn_mapping_to_weights = {}
	# Knn mapping done seperately for each individual
	cell_individuals = np.asarray(adata.obs['ind_cov'])
	unique_individuals = np.unique(cell_individuals)
	cell_names = np.asarray(adata.obs['cell_id'])
	cell_names_to_index = {}
	for i, cell_name in enumerate(cell_names):
		cell_names_to_index[cell_name] = i

	standardized_pca_loadings = standardize_columns(adata.obsm['X_pca'])

	# Loop through individuals
	for individual in unique_individuals:
		print(individual)
		# Get indices of cells corresponding to this individual
		individual_indices = np.where(cell_individuals == individual)[0]
		# Simple error checking
		if 2.0*k > len(individual_indices):
			print('assumption error: not enough cells per individual')
			pdb.set_trace()
		# Filter data to strictly cells from this individual
		individual_expression = adata.X[individual_indices, :]
		individual_pca_loadings = standardized_pca_loadings[individual_indices, :20]
		individual_cell_names = cell_names[individual_indices]
		# Simple error checking
		if len(np.unique(cell_individuals[individual_indices])) != 1:
			print('assumption error')
			pdb.set_trace()
		# compute similarity between cells
		if knn_method.startswith('pearson_full'):
			cell_similarity_matrix = create_cell_similarity_matrix_with_pearson_correlation(individual_expression)
		elif knn_method.startswith('euclidean_pca'):
			cell_similarity_matrix = create_cell_similarity_matrix_with_euclidean_distance(individual_pca_loadings)
		# Loop through cells
		for i, cell_name in enumerate(individual_cell_names):
			cells_similarity_vector = cell_similarity_matrix[i,:]
			top_k_cell_indices = cells_similarity_vector.argsort()[-(k+1):][::-1]
			cell_kernel_width = -cells_similarity_vector[top_k_cell_indices[-1]]
			cell_id_to_kernel_width[cell_name] = cell_kernel_width
	cell_kernel_widths = []
	for cell_name in cell_names:
		cell_kernel_widths.append(cell_id_to_kernel_width[cell_name])
	cell_kernel_widths = np.asarray(cell_kernel_widths)
	median_kernel_width = np.median(cell_kernel_widths)

	# Loop through individuals
	for individual in unique_individuals:
		print(individual)
		# Get indices of cells corresponding to this individual
		individual_indices = np.where(cell_individuals == individual)[0]
		# Simple error checking
		if 2.0*k > len(individual_indices):
			print('assumption error: not enough cells per individual')
			pdb.set_trace()
		# Filter data to strictly cells from this individual
		individual_expression = adata.X[individual_indices, :]
		individual_pca_loadings = standardized_pca_loadings[individual_indices, :20]
		individual_cell_names = cell_names[individual_indices]
		# Simple error checking
		if len(np.unique(cell_individuals[individual_indices])) != 1:
			print('assumption error')
			pdb.set_trace()
		# compute similarity between cells
		if knn_method.startswith('pearson_full'):
			cell_similarity_matrix = create_cell_similarity_matrix_with_pearson_correlation(individual_expression)
		elif knn_method.startswith('euclidean_pca'):
			cell_similarity_matrix = create_cell_similarity_matrix_with_euclidean_distance(individual_pca_loadings)
		# Loop through cells
		for i, cell_name in enumerate(individual_cell_names):
			cells_similarity_vector = cell_similarity_matrix[i,:]
			top_k_cell_indices = cells_similarity_vector.argsort()[-((3*k)+1):][::-1]
			if len(top_k_cell_indices) != (3*k+1):
				print('assumption error')
				pdb.set_trace()
			top_k_cell_names = []
			top_k_indices = []
			top_k_cell_weights = []
			if knn_method == 'euclidean_pca_adaptive_gaussian_kernel':
				cell_kernel_width = cell_id_to_kernel_width[cell_name]
			elif knn_method == 'euclidean_pca_median_gaussian_kernel':
				cell_kernel_width = median_kernel_width
			for cell_index in top_k_cell_indices:
				if cell_index == i:
					continue
				top_k_cell_names.append(individual_cell_names[cell_index])
				index = cell_names_to_index[individual_cell_names[cell_index]]
				top_k_indices.append(index)
				distance = -cells_similarity_vector[cell_index]
				weight = np.exp(-np.square(distance/cell_kernel_width))
				top_k_cell_weights.append(weight)
			top_k_cell_names = np.asarray(top_k_cell_names)
			top_k_indices = np.asarray(top_k_indices)
			top_k_cell_weights = np.asarray(top_k_cell_weights)
			# Error checking
			if len(top_k_indices) != 3*k or len(np.unique(top_k_indices)) != 3*k:
				print('assumption error')
				pdb.set_trace()
			if len(top_k_cell_names) != 3*k or len(np.unique(top_k_cell_names)) != 3*k:
				print('assumption error')
				pdb.set_trace()
			# Add info to dictionary
			if cell_name in knn_mapping_to_indices:
				print('assumption error')
				pdb.set_trace()
			if cell_name in knn_mapping_to_cell_ids:
				print('assumption erorr')
				pdb.set_trace()
			knn_mapping_to_indices[cell_name] = top_k_indices
			knn_mapping_to_cell_ids[cell_name] = top_k_cell_names
			knn_mapping_to_weights[cell_name] = top_k_cell_weights
	if len(cell_names) != len(knn_mapping_to_indices):
		print('aasssumptionerror')
		pdb.set_trace()
	return knn_mapping_to_indices, knn_mapping_to_cell_ids, knn_mapping_to_weights, cell_names, cell_id_to_kernel_width

def generate_knn_boosted_expression_data_wrapper(adata, k, knn_method, raw_knn_boosted_expression_file, standardized_knn_boosted_expression_file, knn_mapping_file, knn_mapping_ct_summary_file, knn_boosted_pca_file, knn_boosted_pca_ve_file, num_pcs):
	if knn_method == 'euclidean_pca' or knn_method == 'pearson_full':
		knn_mapping_to_indices, knn_mapping_to_cell_ids, ordered_cell_ids = create_knn_mapping(adata, k, knn_method)
		print_knn_mapping(knn_mapping_to_indices, knn_mapping_to_cell_ids, ordered_cell_ids, knn_mapping_file)
		print_knn_mapping_cell_type_summary(adata, knn_mapping_file, knn_mapping_ct_summary_file)
		print_knn_mapping_batch_summary(adata, knn_mapping_file, knn_mapping_batch_summary_file)
		# Generate summed knn-boosted counts
		generate_knn_boosted_counts(adata, knn_mapping_to_indices, ordered_cell_ids, raw_knn_boosted_expression_file)
	elif knn_method == 'euclidean_pca_adaptive_gaussian_kernel' or knn_method == 'euclidean_pca_median_gaussian_kernel':
		knn_mapping_to_indices, knn_mapping_to_cell_ids, knn_mapping_to_weights, ordered_cell_ids, cell_id_to_kernel_width = create_distance_weighted_knn_mapping(adata, k, knn_method)
		print_distance_weighted_knn_mapping(knn_mapping_to_indices, knn_mapping_to_cell_ids, knn_mapping_to_weights, ordered_cell_ids, cell_id_to_kernel_width, knn_mapping_file)
		print_distance_weighted_knn_mapping_cell_type_summary(adata, knn_mapping_file, knn_mapping_ct_summary_file)
		# Generate summed knn-boosted counts
		generate_distance_weighted_knn_boosted_counts(adata, knn_mapping_to_indices, ordered_cell_ids, knn_mapping_to_weights, raw_knn_boosted_expression_file)
		
	# Standardize summed pseudobulk counts (samples X genes)
	standardize_pseudobulk_counts(raw_knn_boosted_expression_file, standardized_knn_boosted_expression_file)
	generate_pca_scores_and_variance_explained(standardized_knn_boosted_expression_file, num_pcs, knn_boosted_pca_file, knn_boosted_pca_ve_file)



def convert_from_old_donor_ids_to_new_donor_ids(old_donor_ids, donor_id_mapping):
	new_donor_ids = []
	f = open(donor_id_mapping)
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[1] in mapping:
			if data[0] != mapping[data[1]]:
				print('assumption error')
				pdb.set_trace()
		if len(data) != 2:
			print('assumption errror')
			pdb.set_trace()
		mapping[data[1]] = data[0]
	f.close()
	for old_donor_id in old_donor_ids:
		new_donor_ids.append(mapping[old_donor_id])
	new_donor_ids = np.asarray(new_donor_ids)
	if len(np.unique(old_donor_ids)) != len(np.unique(new_donor_ids)):
		print('assumption error')
		pdb.set_trace()
	if len(new_donor_ids) != len(old_donor_ids):
		print('assumption erororo')
		pdb.set_trace()
	return new_donor_ids

def extract_categorical_feature(raw_feature_cov):
	all_features = []
	for raw_feature in raw_feature_cov:
		all_features.append(raw_feature)
	unique_categories = np.unique(np.asarray(all_features))
	mapping = {}
	for i, category in enumerate(unique_categories):
		mapping[category] = i

	cat_feature = np.zeros((len(raw_feature_cov), len(unique_categories) -1))

	for i, raw_feature in enumerate(raw_feature_cov):
		column_num = mapping[raw_feature]
		if column_num == (len(unique_categories) -1):
			continue
		cat_feature[i, column_num] = 1.0
	return(cat_feature)

def regress_out_batch_effects(adata, expression_file, residual_expression_file, pca_file, pca_ve_file, num_pcs):
	expr = np.loadtxt(expression_file)
	#np.save('temp.npy', expr)
	batch_cov = extract_categorical_feature(adata.obs.batch_cov)
	cov = np.hstack((np.asmatrix(adata.obs['n_counts']).T, batch_cov))
	#expr = np.load('temp.npy')
	reg = LinearRegression().fit(cov, expr)
	print('model fit')
	prediction = reg.predict(cov)
	residual = expr - prediction
	np.savetxt(residual_expression_file, np.round(residual, decimals=5), fmt="%s", delimiter='\t')
	generate_pca_scores_and_variance_explained(residual_expression_file, num_pcs, pca_file, pca_ve_file)

def generate_raw_cluster_pseudobulk_expression(raw_x, ordered_pseudobulk_samples, cluster_assignments, ordered_genes):
	num_samples = len(ordered_pseudobulk_samples)
	num_genes = ordered_genes.shape[0]
	pseudobulk_expr = np.zeros((num_samples, num_genes))

	for pseudobulk_sample_num, pseudobulk_sample_name in enumerate(ordered_pseudobulk_samples):
		# Get cell indices corresponding to this pseudobulk sample
		indices = cluster_assignments == pseudobulk_sample_name
		# Fill in pseudobulk expr
		pseudobulk_expr[pseudobulk_sample_num, :] = np.asarray(np.sum(raw_x[indices,:], axis=0))[0,:]
	return pseudobulk_expr

def get_mode_of_string_list(arr):
	unique,pos = np.unique(arr,return_inverse=True)
	counts = np.bincount(pos)
	maxpos = counts.argmax() 
	mode_value = unique[maxpos]
	return mode_value

# Extract covariates of pseudobulk samples from cell level covariates file
def print_pseudobulk_covariate_file_from_cell_covariates(ordered_pseudobulk_samples, adata_obs, cluster_assignments, pseudobulk_covariate_file):
	# Create mapping from cell type to index position of those cell types
	unique_cell_types = np.unique(adata_obs['cg_cov'])
	num_cell_types = len(unique_cell_types)
	cell_type_to_position_mapping = {}
	for i, cell_type in enumerate(unique_cell_types):
		cell_type_to_position_mapping[cell_type] = i

	# Open output file handle
	t = open(pseudobulk_covariate_file, 'w')
	# Print header
	t.write('pseudobulk_sample\tind_cov\tAge\tSex\tpop_cov\tStatus\tSLE_status\tnum_cells\tcg_cov_mode\tct_cov_mode')
	for cell_type in unique_cell_types:
		t.write('\t' + cell_type + '_fraction')
	t.write('\n')

	count = 0
	for pseudobulk_sample in ordered_pseudobulk_samples:
		t.write(pseudobulk_sample)
		pseudobulk_sample_indices = cluster_assignments == pseudobulk_sample
		# Donor Id
		ind_covs = np.unique(adata_obs['ind_cov'][pseudobulk_sample_indices])
		if len(ind_covs) != 1:
			print('assumption error')
			pdb.set_trace()
		donor_id = ind_covs[0]
		t.write('\t' + donor_id)

		# Average age
		age = np.mean(adata_obs['Age'][pseudobulk_sample_indices].astype(float))
		t.write('\t' + str(age))

		# Sex
		sexs = np.unique(adata_obs['Sex'][pseudobulk_sample_indices])
		if len(sexs) != 1:
			print('assumption eroror')
			pdb.set_trace()
		sex = sexs[0]
		t.write('\t' + sex)

		# Pop cov
		pop_covs = np.unique(adata_obs['pop_cov'][pseudobulk_sample_indices])
		if len(pop_covs) != 1:
			print('assumption erroror')
			pdb.set_trace()
		pop_cov = pop_covs[0]
		t.write('\t' + pop_cov)

		# Status
		statuses = np.unique(adata_obs['Status'][pseudobulk_sample_indices])
		if len(statuses) != 1:
			status = get_mode_of_string_list(np.asarray(adata_obs['Status'][pseudobulk_sample_indices]))
		else:
			status = statuses[0]
		t.write('\t' + status)

		# SLE Status 
		sle_statuses = np.unique(adata_obs['SLE_status'][pseudobulk_sample_indices])
		if len(sle_statuses) != 1:
			print('assumption erorororo')
			pdb.set_trace()
		sle_status = sle_statuses[0]
		t.write('\t' + sle_status)

		# Number of cells in cluster
		num_cells_in_cluster = np.sum(pseudobulk_sample_indices)
		t.write('\t' + str(num_cells_in_cluster))

		# Mode cg_cov
		cg_cov_mode = get_mode_of_string_list(np.asarray(adata_obs['cg_cov'][pseudobulk_sample_indices]))
		t.write('\t' + cg_cov_mode)

		# Mode ct_cov
		ct_cov_mode = get_mode_of_string_list(np.asarray(adata_obs['ct_cov'][pseudobulk_sample_indices]).astype(str))
		t.write('\t' + ct_cov_mode)

		# cg_fraction
		cg_fraction = np.zeros(num_cell_types)
		for cg in np.asarray(adata_obs['cg_cov'][pseudobulk_sample_indices]):
			cg_fraction[cell_type_to_position_mapping[cg]] = cg_fraction[cell_type_to_position_mapping[cg]] + 1
		cg_fraction = cg_fraction/np.sum(cg_fraction)
		t.write('\t' + '\t'.join(cg_fraction.astype(str)) + '\n')
	t.close()

def get_dictionary_list_of_genotyped_individuals(genotyped_individuals_file):
	aa = np.loadtxt(genotyped_individuals_file, dtype=str, delimiter='\t')
	dicti = {}
	for ele in aa:
		dicti[ele] = 1
	return dicti


def generate_cluster_pseudobulk_expression_scran_ign(adata, gene_annotation_file, cluster_assignments, min_depth_threshold, genotyped_individuals_file, output_root):
	# Get dictionary list of genotyped individuals
	geno_indi_dicti = get_dictionary_list_of_genotyped_individuals(genotyped_individuals_file)

	# First generate list of ordered pseudobulk samples
	ordered_pseudobulk_samples_raw = sorted(np.unique(cluster_assignments))
	ordered_pseudobulk_samples = []
	used_indis = {}
	for pseudobulk_sample in ordered_pseudobulk_samples_raw:
		indi = pseudobulk_sample.split(':')[0]
		if indi in geno_indi_dicti:
			ordered_pseudobulk_samples.append(pseudobulk_sample)
			used_indis[indi] = 1
	ordered_pseudobulk_samples = np.asarray(ordered_pseudobulk_samples)
	print(str(len(used_indis)) + ' individuals of ' + str(len(geno_indi_dicti)) + ' genotyped indiviudals used in analysis')

	# Then get ordred list of gene ids
	raw_ordered_genes = np.vstack((adata.raw.var.index, adata.raw.var[adata.raw.var.columns[0]])).T

	raw_protein_coding_gene_indices = extract_protein_coding_known_autosomal_genes(adata.raw.var, gene_annotation_file)


	# Filter out genes with zero counts
	gene_counts = np.asarray(np.sum(adata.raw.X,axis=0))[0,:]
	non_zero_genes = gene_counts != 0.0
	temp_raw = adata.raw.X.toarray()
	temp_raw2 = temp_raw[:, non_zero_genes]
	
	np.savetxt(output_root + 'raw_counts.txt', temp_raw2, fmt="%s", delimiter='\t')
	np.savetxt(output_root + 'raw_counts_subset_debug.txt', temp_raw2[:10000,:], fmt="%s", delimiter='\t')


	print('NEED TO DEAL WITH raw_ordered_genes AND raw_protein_coding_gene_indices')

	'''

	# Generate raw cluster pseudobulk expression
	raw_cluster_pseudobulk_expression = generate_raw_cluster_pseudobulk_expression(adata.raw.X, ordered_pseudobulk_samples, cluster_assignments, raw_ordered_genes)

	# Remove pseudobulk samples with fewer than 10K reads (following Josh's paper heree)
	valid_pseudobulk_sample_indices = np.sum(raw_cluster_pseudobulk_expression, axis=1) >= min_depth_threshold
	filtered_raw_cluster_pseudobulk_expression = raw_cluster_pseudobulk_expression[valid_pseudobulk_sample_indices, :]
	filtered_ordered_pseudobulk_samples = np.asarray(ordered_pseudobulk_samples)[valid_pseudobulk_sample_indices]


	# Remove genes with zeros across all samples
	new_gene_indices = np.sum(raw_cluster_pseudobulk_expression, axis=0) > 0.0
	filtered2_raw_cluster_pseudobulk_expression = filtered_raw_cluster_pseudobulk_expression[:, new_gene_indices]
	ordered_genes = raw_ordered_genes[new_gene_indices, :]
	protein_coding_gene_indices = raw_protein_coding_gene_indices[new_gene_indices]

		
	# TMM normalize expression
	#counts_df = pd.DataFrame(np.transpose(filtered_raw_cluster_pseudobulk_expression))
	#tmm_counts = np.transpose(np.asarray(rnaseqnorm.edgeR_cpm(counts_df, normalized_lib_sizes=True, log=False)))
	#log_tmm_counts = np.transpose(np.asarray(rnaseqnorm.edgeR_cpm(counts_df, normalized_lib_sizes=True, log=True)))
	

	# TMMwsp normalization (needs to be run in R)
	np.savetxt('temp_raw.txt', np.transpose(filtered2_raw_cluster_pseudobulk_expression), delimiter='\t', fmt="%s")
	print(np.transpose(filtered2_raw_cluster_pseudobulk_expression).shape)
	os.system('Rscript tmm_normalization.R temp_raw.txt')
	# Re-load in normalized data
	print('start loading')
	tmm_counts = np.transpose(np.loadtxt('temp_raw_cpm.txt',delimiter='\t'))
	log_tmm_counts = np.transpose(np.loadtxt('temp_raw_log_cpm.txt',delimiter='\t'))

	# Filter to genes that pass our filters
	valid_gene_indices = []
	num_samples = tmm_counts.shape[0]
	min_samples = num_samples*.05
	# Require to pass both tmm filter and raw counts filter
	tmm_pass_filter = np.sum(tmm_counts >= .1, axis=0) >= min_samples
	counts_pass_filter = np.sum(filtered2_raw_cluster_pseudobulk_expression >= 5, axis=0) >= min_samples
	# Test for each gene
	for gene_num in range(len(ordered_genes)):
		if tmm_pass_filter[gene_num] == True and counts_pass_filter[gene_num] == True and protein_coding_gene_indices[gene_num] == True:
			valid_gene_indices.append(True)
		else:
			valid_gene_indices.append(False)
	# Array of length number of genes that determines whether gene passed filter
	valid_gene_indices = np.asarray(valid_gene_indices)

	# Remove genes that don't pass our filters
	filtered_log_tmm_counts = log_tmm_counts[:, valid_gene_indices]
	filtered_ordered_genes = ordered_genes[valid_gene_indices, :]

	# THEN STANDARDIZE
	num_filtered_genes = filtered_log_tmm_counts.shape[1]
	for gene_num in range(num_filtered_genes):
		filtered_log_tmm_counts[:, gene_num] = (filtered_log_tmm_counts[:, gene_num] - np.mean(filtered_log_tmm_counts[:, gene_num]))/np.std(filtered_log_tmm_counts[:, gene_num])

	# Save processed gene expression to output file
	# Pseudobulk expression
	pseudobulk_expression_file = output_root + 'log_tmm_normalized_expression.txt'
	np.savetxt(pseudobulk_expression_file, filtered_log_tmm_counts, fmt="%s", delimiter='\t')
	# Gene names
	gene_names_file = output_root + 'gene_names.txt'
	np.savetxt(gene_names_file, filtered_ordered_genes, fmt="%s", delimiter='\t')
	# Sample names
	sample_names_file = output_root + 'sample_names.txt'
	np.savetxt(sample_names_file, filtered_ordered_pseudobulk_samples, fmt="%s", delimiter='\t')
	# Generate pseudobulk covaraite file
	pseudobulk_covariate_file = output_root + 'sample_covariates.txt'
	print_pseudobulk_covariate_file_from_cell_covariates(filtered_ordered_pseudobulk_samples, adata.obs, cluster_assignments, pseudobulk_covariate_file)

	# Run PCA on pseudobulk data
	num_pcs = 100
	pca_file = output_root + 'pca_scores.txt'
	pca_ve_file = output_root + 'pca_pve.txt'
	generate_pca_scores_and_variance_explained(pseudobulk_expression_file, num_pcs, pca_file, pca_ve_file)

	'''

def generate_cluster_pseudobulk_expression(adata, gene_annotation_file, cluster_assignments, min_depth_threshold, genotyped_individuals_file, output_root):
	# Get dictionary list of genotyped individuals
	geno_indi_dicti = get_dictionary_list_of_genotyped_individuals(genotyped_individuals_file)

	# First generate list of ordered pseudobulk samples
	ordered_pseudobulk_samples_raw = sorted(np.unique(cluster_assignments))
	ordered_pseudobulk_samples = []
	used_indis = {}
	for pseudobulk_sample in ordered_pseudobulk_samples_raw:
		indi = pseudobulk_sample.split(':')[0]
		if indi in geno_indi_dicti:
			ordered_pseudobulk_samples.append(pseudobulk_sample)
			used_indis[indi] = 1
	ordered_pseudobulk_samples = np.asarray(ordered_pseudobulk_samples)
	print(str(len(used_indis)) + ' individuals of ' + str(len(geno_indi_dicti)) + ' genotyped indiviudals used in analysis')

	# Then get ordred list of gene ids
	raw_ordered_genes = np.vstack((adata.raw.var.index, adata.raw.var[adata.raw.var.columns[0]])).T

	raw_protein_coding_gene_indices = extract_protein_coding_known_autosomal_genes(adata.raw.var, gene_annotation_file)


	# Generate raw cluster pseudobulk expression
	raw_cluster_pseudobulk_expression = generate_raw_cluster_pseudobulk_expression(adata.raw.X, ordered_pseudobulk_samples, cluster_assignments, raw_ordered_genes)

	# Remove pseudobulk samples with fewer than 10K reads (following Josh's paper heree)
	valid_pseudobulk_sample_indices = np.sum(raw_cluster_pseudobulk_expression, axis=1) >= min_depth_threshold
	filtered_raw_cluster_pseudobulk_expression = raw_cluster_pseudobulk_expression[valid_pseudobulk_sample_indices, :]
	filtered_ordered_pseudobulk_samples = np.asarray(ordered_pseudobulk_samples)[valid_pseudobulk_sample_indices]


	# Remove genes with zeros across all samples
	new_gene_indices = np.sum(raw_cluster_pseudobulk_expression, axis=0) > 0.0
	filtered2_raw_cluster_pseudobulk_expression = filtered_raw_cluster_pseudobulk_expression[:, new_gene_indices]
	ordered_genes = raw_ordered_genes[new_gene_indices, :]
	protein_coding_gene_indices = raw_protein_coding_gene_indices[new_gene_indices]

	'''	
	# TMM normalize expression
	counts_df = pd.DataFrame(np.transpose(filtered_raw_cluster_pseudobulk_expression))
	tmm_counts = np.transpose(np.asarray(rnaseqnorm.edgeR_cpm(counts_df, normalized_lib_sizes=True, log=False)))
	log_tmm_counts = np.transpose(np.asarray(rnaseqnorm.edgeR_cpm(counts_df, normalized_lib_sizes=True, log=True)))
	'''

	# TMMwsp normalization (needs to be run in R)
	np.savetxt('temp_raw.txt', np.transpose(filtered2_raw_cluster_pseudobulk_expression), delimiter='\t', fmt="%s")
	print(np.transpose(filtered2_raw_cluster_pseudobulk_expression).shape)
	os.system('Rscript tmm_normalization.R temp_raw.txt')
	# Re-load in normalized data
	print('start loading')
	tmm_counts = np.transpose(np.loadtxt('temp_raw_cpm.txt',delimiter='\t'))
	log_tmm_counts = np.transpose(np.loadtxt('temp_raw_log_cpm.txt',delimiter='\t'))

	# Filter to genes that pass our filters
	valid_gene_indices = []
	num_samples = tmm_counts.shape[0]
	min_samples = num_samples*.05
	# Require to pass both tmm filter and raw counts filter
	tmm_pass_filter = np.sum(tmm_counts >= .1, axis=0) >= min_samples
	counts_pass_filter = np.sum(filtered2_raw_cluster_pseudobulk_expression >= 5, axis=0) >= min_samples
	# Test for each gene
	for gene_num in range(len(ordered_genes)):
		if tmm_pass_filter[gene_num] == True and counts_pass_filter[gene_num] == True and protein_coding_gene_indices[gene_num] == True:
			valid_gene_indices.append(True)
		else:
			valid_gene_indices.append(False)
	# Array of length number of genes that determines whether gene passed filter
	valid_gene_indices = np.asarray(valid_gene_indices)

	# Remove genes that don't pass our filters
	filtered_log_tmm_counts = log_tmm_counts[:, valid_gene_indices]
	filtered_ordered_genes = ordered_genes[valid_gene_indices, :]

	# THEN STANDARDIZE
	num_filtered_genes = filtered_log_tmm_counts.shape[1]
	for gene_num in range(num_filtered_genes):
		filtered_log_tmm_counts[:, gene_num] = (filtered_log_tmm_counts[:, gene_num] - np.mean(filtered_log_tmm_counts[:, gene_num]))/np.std(filtered_log_tmm_counts[:, gene_num])

	# Save processed gene expression to output file
	# Pseudobulk expression
	pseudobulk_expression_file = output_root + 'log_tmm_normalized_expression.txt'
	np.savetxt(pseudobulk_expression_file, filtered_log_tmm_counts, fmt="%s", delimiter='\t')
	# Gene names
	gene_names_file = output_root + 'gene_names.txt'
	np.savetxt(gene_names_file, filtered_ordered_genes, fmt="%s", delimiter='\t')
	# Sample names
	sample_names_file = output_root + 'sample_names.txt'
	np.savetxt(sample_names_file, filtered_ordered_pseudobulk_samples, fmt="%s", delimiter='\t')
	# Generate pseudobulk covaraite file
	pseudobulk_covariate_file = output_root + 'sample_covariates.txt'
	print_pseudobulk_covariate_file_from_cell_covariates(filtered_ordered_pseudobulk_samples, adata.obs, cluster_assignments, pseudobulk_covariate_file)

	# Run PCA on pseudobulk data
	num_pcs = 100
	pca_file = output_root + 'pca_scores.txt'
	pca_ve_file = output_root + 'pca_pve.txt'
	generate_pca_scores_and_variance_explained(pseudobulk_expression_file, num_pcs, pca_file, pca_ve_file)

def generate_cluster_pseudobulk_expression_tmm_igp(adata, gene_annotation_file, cluster_assignments, min_depth_threshold, genotyped_individuals_file, output_root):
	# Get dictionary list of genotyped individuals
	geno_indi_dicti = get_dictionary_list_of_genotyped_individuals(genotyped_individuals_file)

	# First generate list of ordered pseudobulk samples
	ordered_pseudobulk_samples_raw = sorted(np.unique(cluster_assignments))
	ordered_pseudobulk_samples = []
	used_indis = {}
	for pseudobulk_sample in ordered_pseudobulk_samples_raw:
		indi = pseudobulk_sample.split(':')[0]
		if indi in geno_indi_dicti:
			ordered_pseudobulk_samples.append(pseudobulk_sample)
			used_indis[indi] = 1
	ordered_pseudobulk_samples = np.asarray(ordered_pseudobulk_samples)
	print(str(len(used_indis)) + ' individuals of ' + str(len(geno_indi_dicti)) + ' genotyped indiviudals used in analysis')

	# Then get ordred list of gene ids
	raw_ordered_genes = np.vstack((adata.raw.var.index, adata.raw.var[adata.raw.var.columns[0]])).T

	raw_protein_coding_gene_indices = extract_protein_coding_known_autosomal_genes(adata.raw.var, gene_annotation_file)


	# Generate raw cluster pseudobulk expression
	raw_cluster_pseudobulk_expression = generate_raw_cluster_pseudobulk_expression(adata.raw.X, ordered_pseudobulk_samples, cluster_assignments, raw_ordered_genes)

	# Remove pseudobulk samples with fewer than 10K reads (following Josh's paper heree)
	valid_pseudobulk_sample_indices = np.sum(raw_cluster_pseudobulk_expression, axis=1) >= min_depth_threshold
	filtered_raw_cluster_pseudobulk_expression = raw_cluster_pseudobulk_expression[valid_pseudobulk_sample_indices, :]
	filtered_ordered_pseudobulk_samples = np.asarray(ordered_pseudobulk_samples)[valid_pseudobulk_sample_indices]


	# Remove genes with zeros across all samples
	new_gene_indices = np.sum(raw_cluster_pseudobulk_expression, axis=0) > 0.0
	filtered2_raw_cluster_pseudobulk_expression = filtered_raw_cluster_pseudobulk_expression[:, new_gene_indices]
	ordered_genes = raw_ordered_genes[new_gene_indices, :]
	protein_coding_gene_indices = raw_protein_coding_gene_indices[new_gene_indices]

	# TMM-IGN normalize expression
	counts_df = pd.DataFrame(np.transpose(filtered2_raw_cluster_pseudobulk_expression))
	temp_out = rnaseqnorm.edgeR_cpm(counts_df, normalized_lib_sizes=True, log=False)
	tmm_counts = np.transpose(np.asarray(temp_out))
	norm_df = rnaseqnorm.inverse_normal_transform(temp_out)
	ign_expression = np.transpose(np.asarray(norm_df))
	#log_tmm_counts = np.transpose(np.asarray(rnaseqnorm.edgeR_cpm(counts_df, normalized_lib_sizes=True, log=True)))


	# TMMwsp normalization (needs to be run in R)
	#np.savetxt('temp_raw.txt', np.transpose(filtered2_raw_cluster_pseudobulk_expression), delimiter='\t', fmt="%s")
	#print(np.transpose(filtered2_raw_cluster_pseudobulk_expression).shape)
	#os.system('Rscript tmm_normalization.R temp_raw.txt')
	# Re-load in normalized data
	#print('start loading')
	#tmm_counts = np.transpose(np.loadtxt('temp_raw_cpm.txt',delimiter='\t'))
	#log_tmm_counts = np.transpose(np.loadtxt('temp_raw_log_cpm.txt',delimiter='\t'))

	# Filter to genes that pass our filters
	valid_gene_indices = []
	num_samples = tmm_counts.shape[0]
	min_samples = num_samples*.2
	# Require to pass both tmm filter and raw counts filter
	tmm_pass_filter = np.sum(tmm_counts >= .1, axis=0) >= min_samples
	counts_pass_filter = np.sum(filtered2_raw_cluster_pseudobulk_expression >= 6, axis=0) >= min_samples
	# Test for each gene
	for gene_num in range(len(ordered_genes)):
		if tmm_pass_filter[gene_num] == True and counts_pass_filter[gene_num] == True and protein_coding_gene_indices[gene_num] == True:
			valid_gene_indices.append(True)
		else:
			valid_gene_indices.append(False)
	# Array of length number of genes that determines whether gene passed filter
	valid_gene_indices = np.asarray(valid_gene_indices)

	# Remove genes that don't pass our filters
	filtered_log_tmm_counts = ign_expression[:, valid_gene_indices]
	filtered_ordered_genes = ordered_genes[valid_gene_indices, :]

	# THEN STANDARDIZE
	num_filtered_genes = filtered_log_tmm_counts.shape[1]
	for gene_num in range(num_filtered_genes):
		filtered_log_tmm_counts[:, gene_num] = (filtered_log_tmm_counts[:, gene_num] - np.mean(filtered_log_tmm_counts[:, gene_num]))/np.std(filtered_log_tmm_counts[:, gene_num])

	# Save processed gene expression to output file
	# Pseudobulk expression
	pseudobulk_expression_file = output_root + 'log_tmm_normalized_expression.txt'
	np.savetxt(pseudobulk_expression_file, filtered_log_tmm_counts, fmt="%s", delimiter='\t')
	# Gene names
	gene_names_file = output_root + 'gene_names.txt'
	np.savetxt(gene_names_file, filtered_ordered_genes, fmt="%s", delimiter='\t')
	# Sample names
	sample_names_file = output_root + 'sample_names.txt'
	np.savetxt(sample_names_file, filtered_ordered_pseudobulk_samples, fmt="%s", delimiter='\t')
	# Generate pseudobulk covaraite file
	pseudobulk_covariate_file = output_root + 'sample_covariates.txt'
	print_pseudobulk_covariate_file_from_cell_covariates(filtered_ordered_pseudobulk_samples, adata.obs, cluster_assignments, pseudobulk_covariate_file)

	# Run PCA on pseudobulk data
	num_pcs = 100
	pca_file = output_root + 'pca_scores.txt'
	pca_ve_file = output_root + 'pca_pve.txt'
	generate_pca_scores_and_variance_explained(pseudobulk_expression_file, num_pcs, pca_file, pca_ve_file)



def print_pseudobulk_clustering_mapping_cell_type_summary(adata, ct_summary_file, cluster_name):
	ordered_cell_types = np.asarray(adata.obs['cg_cov'])
	ordered_cell_ids = np.asarray(adata.obs['cell_id'])
	# create mapping from cell id to cell type
	cell_id_to_cell_type = {}
	for i, cell_id in enumerate(ordered_cell_ids):
		cell_id_to_cell_type[cell_id] = ordered_cell_types[i]
	# Get unique list of cell types
	unique_cell_types = sorted(np.unique(ordered_cell_types))
	num_cell_types = len(unique_cell_types)
	cell_type_to_pos = {}
	for i, cell_type in enumerate(unique_cell_types):
		cell_type_to_pos[cell_type] = i
	# Create dictionary to keep track of cell type knn counts
	mapping = {}
	for cell_type in unique_cell_types:
		mapping[cell_type] = np.zeros(num_cell_types)

	cluster_to_indices_mapping = {}
	unique_clusters = np.unique(adata.obs[cluster_name])
	counter = 0
	for cluster_id in unique_clusters:
		counter = counter + 1
		neighbor_indices = adata.obs[cluster_name] == cluster_id
		cluster_to_indices_mapping[cluster_id] = neighbor_indices


	for i, cell_id in enumerate(ordered_cell_ids):
		cell_type = ordered_cell_types[i]
		cell_cluster_id = adata.obs[cluster_name][i]
		cell_index = cell_type_to_pos[cell_type]
		neighbor_indices = cluster_to_indices_mapping[cell_cluster_id]
		for neighbor_cell_type in ordered_cell_types[neighbor_indices]:
			neighbor_position = cell_type_to_pos[neighbor_cell_type]
			mapping[cell_type][neighbor_position] = mapping[cell_type][neighbor_position] + 1
		# Dont coun't the cell itself (only neighbors)
		mapping[cell_type][cell_index] = mapping[cell_type][cell_index] - 1

	t = open(ct_summary_file, 'w')
	t.write('cell_type\t' + '\t'.join(unique_cell_types) + '\n')
	for cell_type in unique_cell_types:
		t.write(cell_type + '\t')
		frac_ct = mapping[cell_type]/np.sum(mapping[cell_type])
		t.write('\t'.join(frac_ct.astype(str)) + '\n')
	t.close()


#####################
# Command line args
######################
input_h5py_file = sys.argv[1]
processed_expression_dir = sys.argv[2]
gene_annotation_file = sys.argv[3]
genotyped_individuals_file = sys.argv[4]


######################
# Filtering parameters
#######################
np.random.seed(0)
sc.settings.verbosity = 3 
transformation_type = 'log_transform'
regress_out_batch = True
expected_cells_per_pseudobulk_sample = 10


######################
# Load in ScanPy data
#######################
# adata = sc.read_h5ad(input_h5py_file)



##################
# Perform joint cell clustering based on provided clusters
##################
#resolution = 3
#adata.obs['individual_louvain_joint_clusters_' + str(resolution)] = np.char.add(np.char.add(np.asarray(adata.obs['ind_cov']).astype(str), ':'), np.asarray(adata.obs['louvain']).astype(str))


##################
# Perform joint cell clustering
##################
#resolution = 3
#sc.pp.neighbors(adata)
#sc.tl.leiden(adata, resolution=resolution)

#adata.obs['individual_leiden_joint_clusters_' + str(resolution)] = np.char.add(np.char.add(np.asarray(adata.obs['ind_cov']).astype(str), ':'), np.asarray(adata.obs['leiden']).astype(str))

'''
##################
# Perform cell clustering seperately in each individual
##################
unique_individuals = sorted(np.unique(adata.obs['ind_cov']))
num_cells = len(adata.obs['ind_cov'])
#adata.obs['kmeans10'] = np.asarray(['unassigned']*num_cells)
cluster_assignments = np.asarray(['unassigned']*num_cells,dtype='<U40')

resolution = 2.5
# Loop through individuals
for individual in unique_individuals:
	# Get cell indices corresponding to this individual
	cell_indices = adata.obs['ind_cov'] == individual
	# Number of cells in this indiviudal
	num_cells_per_indi = sum(cell_indices)
	# Make anndata object for just this individual
	adata_indi = adata[cell_indices, :]
	# Construct neighborhood graph for cells from this indivudal
	sc.pp.neighbors(adata_indi)
	# Perform leiden clustering
	sc.tl.leiden(adata_indi, resolution=resolution)
	# Get leiden cluster assignemnts
	leiden_cluster_assignments = adata_indi.obs['leiden']
	# Add to global vector of assignments
	cluster_assignments[cell_indices] = np.char.add(individual + ':', leiden_cluster_assignments.astype(str))
	# Delete adata_indi from memory
	del adata_indi
adata.obs['individual_leiden_no_cap_clusters_' + str(resolution)] = cluster_assignments
'''



# Save in temporary adata object
temp_h5_output_file = processed_expression_dir + 'scanpy_temp4.h5ad'
#adata.write(temp_h5_output_file)
adata = sc.read_h5ad(temp_h5_output_file)




#######################
# Add cell id to covariate file
#######################
adata.obs['cell_id'] = adata.obs.index


#######################
# Create cell type summary of clustering file
#######################
resolution = 2.5
clustering_ct_summary_file = processed_expression_dir + 'clustering_leiden_no_cap_resolution_' + str(resolution) + '_cell_type_summary.txt'
# print_pseudobulk_clustering_mapping_cell_type_summary(adata, clustering_ct_summary_file, 'individual_leiden_no_cap_clusters_' + str(resolution))


#######################
# Save Covariate Info
#######################
covariate_output_file = processed_expression_dir + 'cell_covariates.txt'
#np.savetxt(covariate_output_file, adata.obs, fmt="%s", delimiter='\t', header='\t'.join(adata.obs.columns), comments='')
'''


#######################
# Save Expression PCs
#######################
pc_output_file = processed_expression_dir + 'cell_expression_pcs.txt'
np.savetxt(pc_output_file, adata.obsm['X_pca'], fmt="%s", delimiter='\t', comments='')


#######################
# Save Expression PCs variance
#######################
pc_pve_output_file = processed_expression_dir + 'cell_expression_pc_percent_variance_explained.txt'
np.savetxt(pc_pve_output_file, adata.uns['pca']['variance_ratio'], fmt="%s", delimiter='\t', comments='')

#######################
# Save UMAP loadings
#######################
umap_output_file = processed_expression_dir + 'cell_expression_umaps.txt'
np.savetxt(umap_output_file, adata.obsm['X_umap'], fmt="%s", delimiter='\t', comments='')
'''


##################
# Generate cluster-pseudobulk expression
##################
min_depth_threshold = 10000.0
resolution = 2.5
output_root = processed_expression_dir + 'cluster_pseudobulk_leiden_no_cap_' + str(resolution) + '_'
#generate_cluster_pseudobulk_expression(adata, gene_annotation_file, adata.obs['individual_leiden_no_cap_clusters_' + str(resolution)], min_depth_threshold, genotyped_individuals_file, output_root)



##################
# Generate cluster-pseudobulk expression
##################
min_depth_threshold = 50000.0
resolution = 2.5
output_root = processed_expression_dir + 'cluster_tmm_ign_pseudobulk_leiden_no_cap_' + str(resolution) + '_'
# generate_cluster_pseudobulk_expression_tmm_igp(adata, gene_annotation_file, adata.obs['individual_leiden_no_cap_clusters_' + str(resolution)], min_depth_threshold, genotyped_individuals_file, output_root)


##################
# Generate cluster-pseudobulk expression
##################
min_depth_threshold = 50000.0
resolution = 2.5
output_root = processed_expression_dir + 'cluster_scran_ign_pseudobulk_leiden_no_cap_' + str(resolution) + '_'
generate_cluster_pseudobulk_expression_scran_ign(adata, gene_annotation_file, adata.obs['individual_leiden_no_cap_clusters_' + str(resolution)], min_depth_threshold, genotyped_individuals_file, output_root)




































































'''
######################
# Convert from old donor ids to new donors ids
######################
new_donor_ids = convert_from_old_donor_ids_to_new_donor_ids(np.asarray(adata.obs['ind_cov']), genotype_id_mapping_file)
adata.obs['ind_cov'] = new_donor_ids
# Filter out two samples we don't have genotype data for (I know.. a bit hacky)
adata = adata[adata.obs.ind_cov != 'IGTB1540_IGTB1540', :]
adata = adata[adata.obs.ind_cov != '1251_1251', :]

######################
# Filter data
#######################
# Standard filtering
sc.pp.filter_cells(adata, min_genes=300)
sc.pp.filter_genes(adata, min_cells=200)

#######################
# Calculate qc metrics
#######################
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)
# More Standard filtering
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
# Extract which genes are protein-coding, known, autosomal genes
gene_indices = extract_protein_coding_known_autosomal_genes(adata.var, gene_annotation_file)
adata.var['protein_coding_known_autosomal'] = gene_indices
#adata = adata[:, adata.var.protein_coding_known_autosomal==True]

#######################
# Save un-normalized (raw) expression data
#######################
adata.raw = adata

######################
# Normalize data
#######################
if transformation_type == 'log_transform':
	sc.pp.normalize_total(adata, target_sum=1e4)
	sc.pp.log1p(adata)
	sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
	# Filter out highly variable genes
	adata = adata[:, adata.var.highly_variable]
	# If we want to regress things out
	if regress_out_batch == True:
		#sc.pp.regress_out(adata, ['batch_cov'])
		print('start combat')
		sc.pp.combat(adata, key='batch_cov')
		print('end combat')
		print('start regress out total counts and pct counts')
		sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
		print('end regress out total counts and pct counts')
	sc.pp.scale(adata, max_value=10)
elif transformation_type == 'pearson_residual':
	adata = calculate_deviance_residuals(adata)
else:
	print('transformation type: ' + transformation_type + 'currently not implemented')
	pdb.set_trace()


##################
# Run PCA
##################
sc.tl.pca(adata, svd_solver='arpack')


##################
# Run UMAP on full data
##################
sc.pp.neighbors(adata)
sc.tl.umap(adata)

##################
# Perform cell clustering seperately in each individual
##################
unique_individuals = sorted(np.unique(adata.obs['ind_cov']))
num_cells = len(adata.obs['ind_cov'])
#adata.obs['kmeans10'] = np.asarray(['unassigned']*num_cells)
kmeans10_assignments = np.asarray(['unassigned']*num_cells,dtype='<U40')
total_num_clusters = 0

# Loop through individuals
for individual in unique_individuals:
	# Get cell indices corresponding to this individual
	cell_indices = adata.obs['ind_cov'] == individual
	num_cells_per_indi = sum(cell_indices)
	# Get number of clusters in this individual assuming we want a specified number of cells per sample
	num_clusters = int(np.floor(num_cells_per_indi/expected_cells_per_pseudobulk_sample))
	total_num_clusters = total_num_clusters + num_clusters
	# Do KMeans clustering
	kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(adata.obsm['X_pca'][cell_indices,:])
	# Save kmeans clustering assignments for this individual
	kmeans10_assignments[cell_indices] = np.char.add(individual + ':', kmeans.labels_.astype(str))
# Put in adata object
adata.obs['kmeans10'] = kmeans10_assignments
print(total_num_clusters)
print(len(np.unique(kmeans10_assignments)))


'''

'''
##################
# Load in intermediate data
##################
temp_h5_output_file = processed_expression_dir + 'scanpy_temp4.h5ad'
#adata.write(temp_h5_output_file)
adata = sc.read_h5ad(temp_h5_output_file)



##################
# Generate cluster-pseudobulk expression
##################
output_root = processed_expression_dir + 'cluster_pseudobulk_kmeans10_'
generate_cluster_pseudobulk_expression(adata, adata.obs['kmeans10'], output_root)





#######################
# Save Covariate Info
#######################
adata.obs['cell_id'] = adata.obs.index
covariate_output_file = processed_expression_dir + 'cell_covariates.txt'
np.savetxt(covariate_output_file, adata.obs, fmt="%s", delimiter='\t', header='\t'.join(adata.obs.columns), comments='')

#######################
# Save Expression PCs
#######################
pc_output_file = processed_expression_dir + 'cell_expression_pcs.txt'
np.savetxt(pc_output_file, adata.obsm['X_pca'], fmt="%s", delimiter='\t', comments='')


#######################
# Save Expression PCs variance
#######################
pc_pve_output_file = processed_expression_dir + 'cell_expression_pc_percent_variance_explained.txt'
np.savetxt(pc_pve_output_file, adata.uns['pca']['variance_ratio'], fmt="%s", delimiter='\t', comments='')

#######################
# Save UMAP loadings
#######################
umap_output_file = processed_expression_dir + 'cell_expression_umaps.txt'
np.savetxt(umap_output_file, adata.obsm['X_umap'], fmt="%s", delimiter='\t', comments='')

'''

















'''

###########################
# Identify highly variable genes
###########################
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

######################
# Run PCA on data
#######################
sc.tl.pca(adata, svd_solver='arpack')

#######################
# Save Gene IDs
#######################
if random_subset == True:
	gene_id_output_file = processed_expression_dir + 'single_cell_expression_sle_individuals_random_subset_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform' + '_regress_out_batch_' + regress_out_batch_string + '_gene_ids.txt'
else:
	gene_id_output_file = processed_expression_dir + 'single_cell_expression_sle_individuals_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform' + '_regress_out_batch_' + regress_out_batch_string  + '_gene_ids.txt'
np.savetxt(gene_id_output_file, np.vstack((adata.var.index, adata.var[adata.var.columns[0]])).T, fmt="%s", delimiter='\t')

#######################
# Save standardized expression data
#######################
if random_subset == True:
	expression_output_file = processed_expression_dir + 'single_cell_expression_sle_individuals_random_subset_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform'  + '_regress_out_batch_' + regress_out_batch_string + '_standardized.txt'
else:
	expression_output_file = processed_expression_dir + 'single_cell_expression_sle_individuals_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform'  + '_regress_out_batch_' + regress_out_batch_string + '_standardized.txt'
#np.savetxt(expression_output_file, adata.X, fmt="%s", delimiter='\t')

#######################
# Save Covariate Info
#######################
adata.obs['cell_id'] = adata.obs.index
if random_subset == True:
	covariate_output_file = processed_expression_dir + 'cell_covariates_sle_individuals_random_subset_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform' + '_regress_out_batch_' + regress_out_batch_string  + '.txt'
else:
	covariate_output_file = processed_expression_dir + 'cell_covariates_sle_individuals_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform' + '_regress_out_batch_' + regress_out_batch_string + '.txt'
np.savetxt(covariate_output_file, adata.obs, fmt="%s", delimiter='\t', header='\t'.join(adata.obs.columns), comments='')

#######################
# Save PCA Scores
#######################
num_pcs=200
if random_subset == True:
	filtered_cells_pca_file = processed_expression_dir + 'pca_scores_sle_individuals_random_subset_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform' + '_regress_out_batch_' + regress_out_batch_string  + '.txt'
	filtered_cells_pca_ve_file = processed_expression_dir + 'pca_variance_explained_sle_individuals_random_subset_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform' + '_regress_out_batch_' + regress_out_batch_string  + '.txt'
else:
	filtered_cells_pca_file = processed_expression_dir + 'pca_scores_sle_individuals_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform' + '_regress_out_batch_' + regress_out_batch_string  + '.txt'
	filtered_cells_pca_ve_file = processed_expression_dir + 'pca_variance_explained_sle_individuals_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform' + '_regress_out_batch_' + regress_out_batch_string  + '.txt'
generate_pca_scores_and_variance_explained(expression_output_file, num_pcs, filtered_cells_pca_file, filtered_cells_pca_ve_file)

#######################
# Save Output h5 file
#######################
if random_subset == True:
	h5_output_file = processed_expression_dir + 'scanpy_processed_single_cell_data_random_subset_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform' + '_regress_out_batch_' + regress_out_batch_string  + '.h5ad'
else:
	h5_output_file = processed_expression_dir + 'scanpy_processed_single_cell_data_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform' + '_regress_out_batch_' + regress_out_batch_string  + '.h5ad'
adata.write(h5_output_file)

#######################
# Make kNN-boosted expression
#######################
k=30
knn_method = 'euclidean_pca_median_gaussian_kernel'
# Load in data
if random_subset == True:
	processed_single_cell_h5_file = processed_expression_dir + 'scanpy_processed_single_cell_data_random_subset_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'+ '_regress_out_batch_' + regress_out_batch_string  + '.h5ad'
	raw_knn_boosted_expression_file = processed_expression_dir + 'knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string  + '_raw_expression_sle_individuals_random_subset.txt'
	standardized_knn_boosted_expression_file = processed_expression_dir + 'knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_expression_sle_individuals_random_subset_standardized.txt'
	knn_mapping_file = processed_expression_dir + 'knn_mapping_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals_random_subset.txt'
	knn_mapping_ct_summary_file = processed_expression_dir + 'knn_mapping_cell_type_summary_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals_random_subset.txt'
	knn_boosted_pca_file = processed_expression_dir + 'pca_scores_knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals_random_subset.txt'
	knn_boosted_pca_ve_file = processed_expression_dir + 'pca_variance_knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_explained_sle_individuals_random_subset.txt'
else:
	processed_single_cell_h5_file = processed_expression_dir + 'scanpy_processed_single_cell_data_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform'+ '_regress_out_batch_' + regress_out_batch_string + '.h5ad'
	raw_knn_boosted_expression_file = processed_expression_dir + 'knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_raw_expression_sle_individuals.txt'
	standardized_knn_boosted_expression_file = processed_expression_dir + 'knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_expression_sle_individuals_standardized.txt'
	knn_mapping_file = processed_expression_dir + 'knn_mapping_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals.txt'
	knn_mapping_ct_summary_file = processed_expression_dir + 'knn_mapping_cell_type_summary_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals.txt'
	knn_boosted_pca_file = processed_expression_dir + 'pca_scores_knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals.txt'
	knn_boosted_pca_ve_file = processed_expression_dir + 'pca_variance_knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_explained_sle_individuals.txt'
print('start loading')
print(random_subset)
adata = sc.read_h5ad(processed_single_cell_h5_file)
num_pcs=200
print('start')
#generate_knn_boosted_expression_data_wrapper(adata, k, knn_method, raw_knn_boosted_expression_file, standardized_knn_boosted_expression_file, knn_mapping_file, knn_mapping_ct_summary_file, knn_boosted_pca_file, knn_boosted_pca_ve_file, num_pcs)


# Regress out batch
# Currently assumes no random subsetting
standardized_knn_boosted_residaul_expression_file = processed_expression_dir + 'knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_residual_expression_sle_individuals_standardized.txt'
knn_boosted_residual_expression_pca_file = processed_expression_dir + 'pca_scores_knn_boosted_residual_expression_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals.txt'
knn_boosted_residual_expression_pca_ve_file = processed_expression_dir + 'pca_variance_knn_boosted_residual_expression_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_explained_sle_individuals.txt'
num_pcs = 200
#regress_out_batch_effects(adata, standardized_knn_boosted_expression_file, standardized_knn_boosted_residaul_expression_file, knn_boosted_residual_expression_pca_file, knn_boosted_residual_expression_pca_ve_file, num_pcs)

# Regress out batch
# Currently assumes no random subsetting
standardized_residaul_expression_file = processed_expression_dir + 'single_cell_residual_expression_sle_individuals_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform'  + '_regress_out_batch_' + regress_out_batch_string + '_standardized.txt'
residual_expression_pca_file = processed_expression_dir + 'pca_scores_residual_expression_sle_individuals_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform' + '_regress_out_batch_' + regress_out_batch_string  + '.txt'
residual_expression_pca_ve_file = processed_expression_dir + 'pca_variance_residual_expression_explained_sle_individuals_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform' + '_regress_out_batch_' + regress_out_batch_string  + '.txt'
num_pcs = 200
regress_out_batch_effects(adata, expression_output_file, standardized_residaul_expression_file, residual_expression_pca_file, residual_expression_pca_ve_file, num_pcs)



'''





'''
#######################
# Make kNN-boosted expression
#######################
k=30
knn_method = 'euclidean_pca_adaptive_gaussian_kernel'
# Load in data
if random_subset == True:
	processed_single_cell_h5_file = processed_expression_dir + 'scanpy_processed_single_cell_data_random_subset_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'+ '_regress_out_batch_' + regress_out_batch_string  + '.h5ad'
	raw_knn_boosted_expression_file = processed_expression_dir + 'knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string  + '_raw_expression_sle_individuals_random_subset.txt'
	standardized_knn_boosted_expression_file = processed_expression_dir + 'knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_expression_sle_individuals_random_subset_standardized.txt'
	knn_mapping_file = processed_expression_dir + 'knn_mapping_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals_random_subset.txt'
	knn_mapping_ct_summary_file = processed_expression_dir + 'knn_mapping_cell_type_summary_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals_random_subset.txt'
	knn_boosted_pca_file = processed_expression_dir + 'pca_scores_knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals_random_subset.txt'
	knn_boosted_pca_ve_file = processed_expression_dir + 'pca_variance_knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_explained_sle_individuals_random_subset.txt'
else:
	processed_single_cell_h5_file = processed_expression_dir + 'scanpy_processed_single_cell_data_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform'+ '_regress_out_batch_' + regress_out_batch_string + '.h5ad'
	raw_knn_boosted_expression_file = processed_expression_dir + 'knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_raw_expression_sle_individuals.txt'
	standardized_knn_boosted_expression_file = processed_expression_dir + 'knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_expression_sle_individuals_standardized.txt'
	knn_mapping_file = processed_expression_dir + 'knn_mapping_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals.txt'
	knn_mapping_ct_summary_file = processed_expression_dir + 'knn_mapping_cell_type_summary_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals.txt'
	knn_boosted_pca_file = processed_expression_dir + 'pca_scores_knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals.txt'
	knn_boosted_pca_ve_file = processed_expression_dir + 'pca_variance_knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_explained_sle_individuals.txt'
print('start loading')
print(random_subset)
adata = sc.read_h5ad(processed_single_cell_h5_file)
num_pcs=200
print('start')
generate_knn_boosted_expression_data_wrapper(adata, k, knn_method, raw_knn_boosted_expression_file, standardized_knn_boosted_expression_file, knn_mapping_file, knn_mapping_ct_summary_file, knn_boosted_pca_file, knn_boosted_pca_ve_file, num_pcs)

#######################
# Make kNN-boosted expression
#######################
k=90
knn_method = 'euclidean_pca'

# Load in data
if random_subset == True:
	processed_single_cell_h5_file = processed_expression_dir + 'scanpy_processed_single_cell_data_random_subset_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'+ '_regress_out_batch_' + regress_out_batch_string  + '.h5ad'
	raw_knn_boosted_expression_file = processed_expression_dir + 'knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string  + '_raw_expression_sle_individuals_random_subset.txt'
	standardized_knn_boosted_expression_file = processed_expression_dir + 'knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_expression_sle_individuals_random_subset_standardized.txt'
	knn_mapping_file = processed_expression_dir + 'knn_mapping_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals_random_subset.txt'
	knn_mapping_ct_summary_file = processed_expression_dir + 'knn_mapping_cell_type_summary_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals_random_subset.txt'
	knn_boosted_pca_file = processed_expression_dir + 'pca_scores_knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals_random_subset.txt'
	knn_boosted_pca_ve_file = processed_expression_dir + 'pca_variance_knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_explained_sle_individuals_random_subset.txt'
else:
	processed_single_cell_h5_file = processed_expression_dir + 'scanpy_processed_single_cell_data_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform'+ '_regress_out_batch_' + regress_out_batch_string + '.h5ad'
	raw_knn_boosted_expression_file = processed_expression_dir + 'knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_raw_expression_sle_individuals.txt'
	standardized_knn_boosted_expression_file = processed_expression_dir + 'knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_expression_sle_individuals_standardized.txt'
	knn_mapping_file = processed_expression_dir + 'knn_mapping_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals.txt'
	knn_mapping_ct_summary_file = processed_expression_dir + 'knn_mapping_cell_type_summary_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals.txt'
	knn_boosted_pca_file = processed_expression_dir + 'pca_scores_knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_sle_individuals.txt'
	knn_boosted_pca_ve_file = processed_expression_dir + 'pca_variance_knn_boosted_k_' + str(k) + '_' + knn_method + '_regress_out_batch_' + regress_out_batch_string + '_explained_sle_individuals.txt'
print('start loading')
print(random_subset)
adata = sc.read_h5ad(processed_single_cell_h5_file)
num_pcs=200
print('start')
generate_knn_boosted_expression_data_wrapper(adata, k, knn_method, raw_knn_boosted_expression_file, standardized_knn_boosted_expression_file, knn_mapping_file, knn_mapping_ct_summary_file, knn_boosted_pca_file, knn_boosted_pca_ve_file, num_pcs)

'''




































































############################
# OLD (Not currently used)
############################










'''
#######################
# Make Pseudo-bulk expression data (seperate sample for each cell-type, individual)
#######################
if random_subset == True:
	processed_single_cell_h5_file = processed_expression_dir + 'scanpy_processed_single_cell_data_random_subset.h5ad'
	raw_pseudobulk_expression_file = processed_expression_dir + 'pseudobulk_raw_expression_sle_individuals_random_subset.txt'
	standardized_pseudobulk_expression_file = processed_expression_dir + 'pseudobulk_expression_sle_individuals_random_subset_standardized.txt'
	pseudobulk_covariate_file = processed_expression_dir + 'pseudobulk_covariates_sle_individuals_random_subset.txt'
else:
	processed_single_cell_h5_file = processed_expression_dir + 'scanpy_processed_single_cell_data.h5ad'
	raw_pseudobulk_expression_file = processed_expression_dir + 'pseudobulk_raw_expression_sle_individuals.txt'
	standardized_pseudobulk_expression_file = processed_expression_dir + 'pseudobulk_expression_sle_individuals_standardized.txt'
	pseudobulk_covariate_file = processed_expression_dir + 'pseudobulk_covariates_sle_individuals.txt'
adata = sc.read_h5ad(processed_single_cell_h5_file)
generate_pseudobulk_expression_data_wrapper(adata, raw_pseudobulk_expression_file, standardized_pseudobulk_expression_file, pseudobulk_covariate_file)
#######################
# Save Pseudobulk PCA Scores
#######################
num_pcs=200
if random_subset == True:
	filtered_cells_pca_file = processed_expression_dir + 'pca_scores_pseudobulk_sle_individuals_random_subset.txt'
	filtered_cells_pca_ve_file = processed_expression_dir + 'pca_variance_pseudobulk_explained_sle_individuals_random_subset.txt'
else:
	filtered_cells_pca_file = processed_expression_dir + 'pca_scores_pseudobulk_sle_individuals.txt'
	filtered_cells_pca_ve_file = processed_expression_dir + 'pca_variance_pseudobulk_explained_sle_individuals.txt'
generate_pca_scores_and_variance_explained(standardized_pseudobulk_expression_file, num_pcs, filtered_cells_pca_file, filtered_cells_pca_ve_file)
'''







'''
#######################
# Make Pseudo-bulk expression data (seperate sample for individual)
#######################
if random_subset == True:
	processed_single_cell_h5_file = processed_expression_dir + 'scanpy_processed_single_cell_data_random_subset_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform'+ '_regress_out_batch_' + regress_out_batch_string + '.h5ad'
	raw_pseudobulk_expression_file = processed_expression_dir + 'pseudobulk_per_individual_raw_expression_sle_individuals_random_subset.txt'
	standardized_pseudobulk_expression_file = processed_expression_dir + 'pseudobulk_per_individual_expression_sle_individuals_random_subset_standardized.txt'
	pseudobulk_covariate_file = processed_expression_dir + 'pseudobulk_per_individual_covariates_sle_individuals_random_subset.txt'
else:
	processed_single_cell_h5_file = processed_expression_dir  + 'scanpy_processed_single_cell_data_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform'+ '_regress_out_batch_' + regress_out_batch_string + '.h5ad'
	raw_pseudobulk_expression_file = processed_expression_dir + 'pseudobulk_per_individual_raw_expression_sle_individuals.txt'
	standardized_pseudobulk_expression_file = processed_expression_dir + 'pseudobulk_per_individual_expression_sle_individuals_standardized.txt'
	pseudobulk_covariate_file = processed_expression_dir + 'pseudobulk_per_individual_covariates_sle_individuals.txt'
adata = sc.read_h5ad(processed_single_cell_h5_file)
generate_pseudobulk_per_individual_expression_data_wrapper(adata, raw_pseudobulk_expression_file, standardized_pseudobulk_expression_file, pseudobulk_covariate_file)
#######################
# Save Pseudobulk PCA Scores
#######################
num_pcs=200
if random_subset == True:
	filtered_cells_pca_file = processed_expression_dir + 'pca_scores_pseudobulk_per_individual_sle_individuals_random_subset.txt'
	filtered_cells_pca_ve_file = processed_expression_dir + 'pca_variance_pseudobulk_per_individual_explained_sle_individuals_random_subset.txt'
else:
	filtered_cells_pca_file = processed_expression_dir + 'pca_scores_pseudobulk_per_individual_sle_individuals.txt'
	filtered_cells_pca_ve_file = processed_expression_dir + 'pca_variance_pseudobulk_per_individual_explained_sle_individuals.txt'
generate_pca_scores_and_variance_explained(standardized_pseudobulk_expression_file, num_pcs, filtered_cells_pca_file, filtered_cells_pca_ve_file)
'''

'''
#######################
# Make cell type specific expression data
#######################
# Load in data
if random_subset == True:
	processed_single_cell_h5_file = processed_expression_dir + 'scanpy_processed_single_cell_data_random_subset_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'  + '.h5ad'
else:
	processed_single_cell_h5_file = processed_expression_dir + 'scanpy_processed_single_cell_data_min_expressed_cells_' + str(min_fraction_of_cells) + '_' + transformation_type + '_transform'  + '.h5ad'
adata = sc.read_h5ad(processed_single_cell_h5_file)
# Get unique cell types
unique_cell_types = np.unique(adata.obs.ct_cov)
np.savetxt(processed_expression_dir + 'cell_types.txt', unique_cell_types, fmt="%s")
# Loop through cell types
for cell_type in unique_cell_types:
	print(cell_type)
	printible_cell_type = '_'.join(cell_type.split(' '))
	if random_subset == False:
		standardized_cell_type_expression_file = processed_expression_dir + printible_cell_type + '_single_cell_expression_sle_individuals_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'  + '_standardized.txt'
		raw_cell_type_expression_file = processed_expression_dir + printible_cell_type + '_single_cell_expression_sle_individuals_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'  + '_raw.txt'
		centered_raw_cell_type_expression_file = processed_expression_dir + printible_cell_type + '_single_cell_expression_sle_individuals_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'  + '_raw_centered.txt'
		standardized_cell_type_covariate_file = processed_expression_dir + printible_cell_type + '_cell_covariates_sle_individuals_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'  + '.txt'
		cell_type_pca_file = processed_expression_dir + printible_cell_type + '_pca_scores_sle_individuals_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'  + '.txt'
		cell_type_pca_pve_file = processed_expression_dir + printible_cell_type + '_pca_variance_explained_sle_individuals_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'  + '.txt'
		cell_type_raw_pca_file = processed_expression_dir + printible_cell_type + '_raw_pca_scores_sle_individuals_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'  + '.txt'
		cell_type_raw_pca_pve_file = processed_expression_dir + printible_cell_type + '_raw_pca_variance_explained_sle_individuals_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'  + '.txt'
	else:
		standardized_cell_type_expression_file = processed_expression_dir + printible_cell_type + '_single_cell_expression_sle_individuals_random_subset_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'  + '_standardized.txt'
		raw_cell_type_expression_file = processed_expression_dir + printible_cell_type + '_single_cell_expression_sle_individuals_random_subset_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'  + '_raw.txt'
		centered_raw_cell_type_expression_file = processed_expression_dir + printible_cell_type + '_single_cell_expression_sle_individuals_random_subset_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'  + '_raw_centered.txt'
		standardized_cell_type_covariate_file = processed_expression_dir + printible_cell_type + '_cell_covariates_sle_individuals_random_subset_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'  + '.txt'
		cell_type_pca_file = processed_expression_dir + printible_cell_type + '_pca_scores_sle_individuals_random_subset_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'  + '.txt'
		cell_type_pca_pve_file = processed_expression_dir + printible_cell_type + '_pca_variance_explained_sle_individuals_random_subset_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'  + '.txt'
		cell_type_raw_pca_file = processed_expression_dir + printible_cell_type + '_raw_pca_scores_sle_individuals_random_subset_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'  + '.txt'
		cell_type_raw_pca_pve_file = processed_expression_dir + printible_cell_type + '_raw_pca_variance_explained_sle_individuals_random_subset_min_expressed_cells_' + str(min_fraction_of_cells)+ '_' + transformation_type + '_transform'  + '.txt'
	generate_cell_type_sc_expression_data_wrapper(cell_type, adata, standardized_cell_type_expression_file, standardized_cell_type_covariate_file, raw_cell_type_expression_file, centered_raw_cell_type_expression_file)
	num_pcs=200
	generate_pca_scores_and_variance_explained(standardized_cell_type_expression_file, num_pcs, cell_type_pca_file, cell_type_pca_pve_file)
	generate_pca_scores_and_variance_explained(centered_raw_cell_type_expression_file, num_pcs, cell_type_raw_pca_file, cell_type_raw_pca_pve_file)
'''

'''
#######################
# Make Pseudo-bulk expression data
#######################
if random_subset == True:
	processed_single_cell_h5_file = processed_expression_dir + 'scanpy_processed_single_cell_data_random_subset.h5ad'
	raw_pseudobulk_expression_file = processed_expression_dir + 'pseudobulk_raw_expression_sle_individuals_random_subset.txt'
	standardized_pseudobulk_expression_file = processed_expression_dir + 'pseudobulk_expression_sle_individuals_random_subset_standardized.txt'
	pseudobulk_covariate_file = processed_expression_dir + 'pseudobulk_covariates_sle_individuals_random_subset.txt'
else:
	processed_single_cell_h5_file = processed_expression_dir + 'scanpy_processed_single_cell_data.h5ad'
	raw_pseudobulk_expression_file = processed_expression_dir + 'pseudobulk_raw_expression_sle_individuals.txt'
	standardized_pseudobulk_expression_file = processed_expression_dir + 'pseudobulk_expression_sle_individuals_standardized.txt'
	pseudobulk_covariate_file = processed_expression_dir + 'pseudobulk_covariates_sle_individuals.txt'
adata = sc.read_h5ad(processed_single_cell_h5_file)
generate_pseudobulk_expression_data_wrapper(adata, raw_pseudobulk_expression_file, standardized_pseudobulk_expression_file, pseudobulk_covariate_file)
#######################
# Save Pseudobulk PCA Scores
#######################
num_pcs=200
if random_subset == True:
	filtered_cells_pca_file = processed_expression_dir + 'pca_scores_pseudobulk_sle_individuals_random_subset.txt'
	filtered_cells_pca_ve_file = processed_expression_dir + 'pca_variance_pseudobulk_explained_sle_individuals_random_subset.txt'
else:
	filtered_cells_pca_file = processed_expression_dir + 'pca_scores_pseudobulk_sle_individuals.txt'
	filtered_cells_pca_ve_file = processed_expression_dir + 'pca_variance_pseudobulk_explained_sle_individuals.txt'
generate_pca_scores_and_variance_explained(standardized_pseudobulk_expression_file, num_pcs, filtered_cells_pca_file, filtered_cells_pca_ve_file)
'''
