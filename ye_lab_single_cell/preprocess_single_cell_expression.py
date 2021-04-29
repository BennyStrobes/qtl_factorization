import numpy as np 
import os
import sys
import pdb
import h5py
from sklearn.decomposition import PCA
import scanpy as sc
from anndata import AnnData
from sklearn.metrics.pairwise import euclidean_distances



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
	uuu, sss, vh = np.linalg.svd(np.transpose(X), full_matrices=False)
	svd_loadings = np.transpose(vh)[:,:num_pcs]

	# Save to output file
	np.savetxt(filtered_cells_pca_file, svd_loadings, fmt="%s", delimiter='\t')

	# Compute variance explained
	ve = (np.square(sss)/np.sum(np.square(sss)))[:num_pcs]
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

	np.savetxt(standardized_pseudobulk_expression_file, standardized_pseudobulk_expression, fmt="%s", delimiter='\t')


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
		#print_knn_mapping_batch_summary(adata, knn_mapping_file, knn_mapping_batch_summary_file)
		'''
		# Temp (for speed of writing code)
		#knn_mapping_to_indices, knn_mapping_to_cell_ids, ordered_cell_ids = temp_load_knn_mapping(knn_mapping_file)
		'''
		# Generate summed knn-boosted counts
		#generate_knn_boosted_counts(adata, knn_mapping_to_indices, ordered_cell_ids, raw_knn_boosted_expression_file)
	elif knn_method == 'euclidean_pca_adaptive_gaussian_kernel' or knn_method == 'euclidean_pca_median_gaussian_kernel':
		knn_mapping_to_indices, knn_mapping_to_cell_ids, knn_mapping_to_weights, ordered_cell_ids, cell_id_to_kernel_width = create_distance_weighted_knn_mapping(adata, k, knn_method)
		print_distance_weighted_knn_mapping(knn_mapping_to_indices, knn_mapping_to_cell_ids, knn_mapping_to_weights, ordered_cell_ids, cell_id_to_kernel_width, knn_mapping_file)

		'''
		# Temp (for speed of writing code)
		knn_mapping_to_indices, knn_mapping_to_cell_ids, knn_mapping_to_weights, ordered_cell_ids, cell_id_to_kernel_width = temp_load_distance_weighted_knn_mapping(knn_mapping_file)
		'''
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

#####################
# Command line args
######################
input_h5py_file = sys.argv[1]
processed_expression_dir = sys.argv[2]
gene_annotation_file = sys.argv[3]
min_fraction_of_cells = float(sys.argv[4])
min_genes = int(sys.argv[5])
transformation_type = sys.argv[6]  # Either deviance_residual or log_transform
genotype_id_mapping_file = sys.argv[7]


######################
# Filtering parameters
#######################
#min_genes = 400
# Min fraction of expressed cells for a gene
#min_fraction_of_cells = .1
# Random subset
random_subset = False
regress_out_batch = False
if regress_out_batch == True:
	regress_out_batch_string = 'True'
else:
	regress_out_batch_string = 'False'
np.random.seed(0)

######################
# Load in ScanPy data
#######################
adata = sc.read_h5ad(input_h5py_file)

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
# Allow option to randomly generate subset of the data
if random_subset == True:
	num_cells = adata.X.shape[0]
	adata.obs['random_subset'] = np.random.uniform(size=num_cells) < (1/5)
	adata = adata[adata.obs.random_subset == True, :]
# Standard filtering
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
# Extract percent mito-counts
mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata.obs['n_counts'] = adata.X.sum(axis=1).A1
# Standard filtering
adata = adata[adata.obs.n_genes < 2500, :]
adata = adata[adata.obs.percent_mito < 0.05, :]
# Limit to protein-coding, known, autosomal genes
gene_indices = extract_protein_coding_known_autosomal_genes(adata.var, gene_annotation_file)
adata.var['protein_coding_known_autosomal'] = gene_indices
adata = adata[:, adata.var.protein_coding_known_autosomal==True]
# Filter rows and columns
print(adata.X.shape)
sc.pp.filter_cells(adata, min_genes=min_genes)
print(adata.X.shape)
sc.pp.filter_genes(adata, min_cells=(adata.X.shape[0])*min_fraction_of_cells)
print(adata.X.shape)

#######################
# Save un-normalized (raw) expression data
#######################
if random_subset == True:
	expression_output_file = processed_expression_dir + 'single_cell_raw_expression_sle_individuals_random_subset_min_expressed_cells_' + str(min_fraction_of_cells) + '.txt'
else:
	expression_output_file = processed_expression_dir + 'single_cell_raw_expression_sle_individuals_min_expressed_cells_' + str(min_fraction_of_cells) + '.txt'
np.savetxt(expression_output_file, adata.X.toarray(), fmt="%s", delimiter='\t')
adata.raw = adata

######################
# Normalize data
#######################
if transformation_type == 'log_transform':
	sc.pp.normalize_total(adata, target_sum=1e4)
	sc.pp.log1p(adata)
	if regress_out_batch == True:
		#sc.pp.regress_out(adata, ['batch_cov'])
		print('start combat')
		sc.pp.combat(adata, key='batch_cov')
		print('end combat')
	sc.pp.scale(adata, max_value=10)
elif transformation_type == 'pearson_residual':
	adata = calculate_deviance_residuals(adata)
else:
	print('transformation type: ' + transformation_type + 'currently not implemented')
	pdb.set_trace()

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
np.savetxt(expression_output_file, adata.X, fmt="%s", delimiter='\t')

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
generate_knn_boosted_expression_data_wrapper(adata, k, knn_method, raw_knn_boosted_expression_file, standardized_knn_boosted_expression_file, knn_mapping_file, knn_mapping_ct_summary_file, knn_boosted_pca_file, knn_boosted_pca_ve_file, num_pcs)




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
