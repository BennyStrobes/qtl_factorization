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



def perform_leiden_clustering_in_each_individual(adata, cluster_resolution):
	unique_individuals = sorted(np.unique(adata.obs['ind_cov']))
	num_cells = len(adata.obs['ind_cov'])
	# Initialize vector to keep track of cluster assignments
	cluster_assignments = np.asarray(['unassigned']*num_cells,dtype='<U40')

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
		sc.tl.leiden(adata_indi, resolution=cluster_resolution)
		# Get leiden cluster assignemnts
		leiden_cluster_assignments = adata_indi.obs['leiden']
		# Add to global vector of assignments
		cluster_assignments[cell_indices] = np.char.add(individual + ':', leiden_cluster_assignments.astype(str))
		# Delete adata_indi from memory
		del adata_indi
	adata.obs['individual_leiden_clusters_' + str(cluster_resolution)] = cluster_assignments
	return adata

def get_dictionary_list_of_genotyped_individuals(genotyped_individuals_file):
	aa = np.loadtxt(genotyped_individuals_file, dtype=str, delimiter='\t')
	dicti = {}
	for ele in aa:
		dicti[ele] = 1
	return dicti

def extract_ordered_list_of_pseudobulk_samples(adata, genotyped_individuals_file, cluster_assignments, min_cells_per_indi):
	#Get dictionary list of genotyped individuals
	geno_indi_dicti = get_dictionary_list_of_genotyped_individuals(genotyped_individuals_file)

	# Compute number of cells per individual
	cells_per_indi = {}
	for indi in np.asarray(np.unique(adata.obs['ind_cov'])):
		cells_per_indi[indi] = sum(adata.obs['ind_cov'] == indi)

	# Generate list of ordered pseudobulk samples
	ordered_pseudobulk_samples_raw = sorted(np.unique(cluster_assignments))
	ordered_pseudobulk_samples = []
	used_indis = {}
	# Now make sure each of these pseudobulk samples is genotyped and comes from an indiviudal with at least min_cells_per_indi
	for pseudobulk_sample in ordered_pseudobulk_samples_raw:
		indi = pseudobulk_sample.split(':')[0]
		if indi in geno_indi_dicti and cells_per_indi[indi] > min_cells_per_indi:
			ordered_pseudobulk_samples.append(pseudobulk_sample)
			used_indis[indi] = 1
	ordered_pseudobulk_samples = np.asarray(ordered_pseudobulk_samples)
	print(str(len(used_indis)) + ' individuals of ' + str(len(geno_indi_dicti)) + ' genotyped indiviudals used in analysis')
	return ordered_pseudobulk_samples


def generate_mean_cluster_pseudobulk_expression(processed_x, ordered_pseudobulk_samples, cluster_assignments, ordered_genes):
	num_samples = len(ordered_pseudobulk_samples)
	num_genes = ordered_genes.shape[0]
	pseudobulk_expr = np.zeros((num_samples, num_genes))

	for pseudobulk_sample_num, pseudobulk_sample_name in enumerate(ordered_pseudobulk_samples):
		# Get cell indices corresponding to this pseudobulk sample
		indices = cluster_assignments == pseudobulk_sample_name
		# Fill in pseudobulk expr
		pseudobulk_expr[pseudobulk_sample_num, :] = np.asarray(np.mean(processed_x[indices,:], axis=0))
	return pseudobulk_expr

def get_mode_of_string_list(arr):
	unique,pos = np.unique(arr,return_inverse=True)
	counts = np.bincount(pos)
	maxpos = counts.argmax() 
	mode_value = unique[maxpos]
	return mode_value

# Extract covariates of pseudobulk samples from cell level covariates file
def print_pseudobulk_covariate_file_from_cell_covariates(ordered_pseudobulk_samples, adata_obs, cluster_assignments, pseudobulk_covariate_file, donor_to_isg_score):
	# Create mapping from cell type to index position of those cell types
	unique_cell_types = np.unique(adata_obs['cg_cov'])
	num_cell_types = len(unique_cell_types)
	cell_type_to_position_mapping = {}
	for i, cell_type in enumerate(unique_cell_types):
		cell_type_to_position_mapping[cell_type] = i

	# Open output file handle
	t = open(pseudobulk_covariate_file, 'w')
	# Print header
	t.write('pseudobulk_sample\tind_cov\tAge\tSex\tpop_cov\tStatus\tSLE_status\tnum_cells\tbatch_cov\tcg_cov_mode\tct_cov_mode')
	for cell_type in unique_cell_types:
		t.write('\t' + cell_type + '_fraction')
	t.write('\tavg_n_genes_by_counts\tavg_log1p_n_genes_by_counts\ttotal_counts\tlog_total_counts\tavg_total_counts\tavg_log_total_counts\tavg_pct_counts_in_top_50_genes\tavg_pct_counts_in_top_100_genes\tavg_pct_counts_in_top_200_genes\tdonor_isg_score')
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

		# Batch cov
		batches = np.unique(adata_obs['batch_cov'][pseudobulk_sample_indices])
		if len(batches) != 1:
			batch = get_mode_of_string_list(np.asarray(adata_obs['batch_cov'][pseudobulk_sample_indices]))
		else:
			batch = batches[0]
		t.write('\t' + batch)


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
		t.write('\t' + '\t'.join(cg_fraction.astype(str)))

		# n_genes_by_counts
		avg_n_genes_by_counts = np.mean(adata.obs['n_genes_by_counts'][pseudobulk_sample_indices])
		t.write('\t' + str(avg_n_genes_by_counts))

		# avg log1p_n_genes_by_counts
		avg_log1p_n_genes_by_counts = np.mean(adata.obs['log1p_n_genes_by_counts'][pseudobulk_sample_indices])
		t.write('\t' + str(avg_log1p_n_genes_by_counts))

		# total_counts
		total_counts = np.sum(adata.obs['total_counts'][pseudobulk_sample_indices])
		t.write('\t' + str(total_counts))	

		# log total_counts
		total_counts = np.sum(adata.obs['total_counts'][pseudobulk_sample_indices])
		t.write('\t' + str(np.log(total_counts)))

		# avg total_counts
		avg_total_counts = np.mean(adata.obs['total_counts'][pseudobulk_sample_indices])
		t.write('\t' + str(avg_total_counts))	

		# avg log total_counts
		avg_log_total_counts = np.mean(adata.obs['log1p_total_counts'][pseudobulk_sample_indices])
		t.write('\t' + str(avg_log_total_counts))	

		# avg pct_counts in top 50 genes
		avg_pct_counts_in_top_50_genes = np.mean(adata.obs['pct_counts_in_top_50_genes'][pseudobulk_sample_indices])
		t.write('\t' + str(avg_pct_counts_in_top_50_genes))	

		# avg pct_counts in top 100 genes
		avg_pct_counts_in_top_100_genes = np.mean(adata.obs['pct_counts_in_top_100_genes'][pseudobulk_sample_indices])
		t.write('\t' + str(avg_pct_counts_in_top_100_genes))	

		# avg pct_counts in top 200 genes
		avg_pct_counts_in_top_200_genes = np.mean(adata.obs['pct_counts_in_top_200_genes'][pseudobulk_sample_indices])
		t.write('\t' + str(avg_pct_counts_in_top_200_genes))

		# Donor isg score
		t.write('\t' + str(donor_to_isg_score[donor_id]) + '\n')
	t.close()


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


def normalize_expression_and_generate_expression_pcs(raw_pseudobulk_expression, sample_level_normalization, gene_level_normalization, num_pcs, pb_expression_output_root):
	# Initialize output normalized expression matrix
	normalized_expression = np.zeros(raw_pseudobulk_expression.shape)

	##################################
	# Perform sample level normalization
	##################################
	if sample_level_normalization == 'qn':
		df = pd.DataFrame(np.transpose(raw_pseudobulk_expression))
		temp_out = rnaseqnorm.normalize_quantiles(df)
		raw_pseudobulk_expression = np.transpose(np.asarray(temp_out))

	##################################
	# Perform gene level normalization
	##################################
	if gene_level_normalization == 'zscore':
		for gene_num in range(normalized_expression.shape[1]):
			normalized_expression[:,gene_num] = (raw_pseudobulk_expression[:, gene_num] - np.mean(raw_pseudobulk_expression[:, gene_num]))/np.std(raw_pseudobulk_expression[:, gene_num])
	elif gene_level_normalization == 'ign':
		# Code from GTEx v8
		# Project each gene onto a gaussian
		df = pd.DataFrame(np.transpose(raw_pseudobulk_expression))
		norm_df = rnaseqnorm.inverse_normal_transform(df)
		normalized_expression = np.transpose(np.asarray(norm_df))
	else:
		print(gene_level_normalization + ' gene level normalization method currently not implemented')
		pdb.set_trace()

	# Save normalized pseudobulk gene expression to output file
	pseudobulk_expression_file = pb_expression_output_root + 'normalized_expression.txt'
	np.savetxt(pseudobulk_expression_file, normalized_expression, fmt="%s", delimiter='\t')

	# Run PCA on pseudobulk data
	pca_file = pb_expression_output_root + 'pca_scores.txt'
	pca_ve_file = pb_expression_output_root + 'pca_pve.txt'
	generate_pca_scores_and_variance_explained(pseudobulk_expression_file, num_pcs, pca_file, pca_ve_file)

def create_donor_to_isg_score_mapping(isg_score_file):
	f = open(isg_score_file)
	mapping = {}
	for line in f:
		line = line.rstrip()
		data = line.split(',')
		donor_id = data[0]
		score = float(data[1])
		mapping[donor_id] = score
	f.close()
	return mapping

#####################
# Command line args
#####################
processed_expression_dir = sys.argv[1]  # Input dir
processed_pseudobulk_expression_dir = sys.argv[2]  # Output dir
genotyped_individuals_file = sys.argv[3]  # File containing list of which individuals are genotyped
cluster_resolution = float(sys.argv[4])
regress_out_batch = sys.argv[5]  # Hyperparameter
isg_score_file = sys.argv[6]

# Load in processed-SC Ann-Data file
input_h5py_file = processed_expression_dir + 'scran_normalization_regress_batch_' + regress_out_batch + '_2.h5ad'
#adata = sc.read_h5ad(input_h5py_file)

#################
# Create mapping from donor id to isg score
################
donor_to_isg_score = create_donor_to_isg_score_mapping(isg_score_file)


##################
# Perform cell clustering seperately in each individual
##################
#adata = perform_leiden_clustering_in_each_individual(adata, cluster_resolution)

# Save in temporary adata object
temp_h5_output_file = processed_pseudobulk_expression_dir + 'scran_normalization_regress_batch_' + regress_out_batch + '_with_individual_leiden_clusters_' + str(cluster_resolution) + '.h5ad'
#adata.write(temp_h5_output_file)
adata = sc.read_h5ad(temp_h5_output_file)


#######################
# Add cell id to covariate file
#######################
adata.obs['cell_id'] = adata.obs.index


#######################
# Create cell type summary output file of clustering file
#######################
clustering_ct_summary_file = processed_pseudobulk_expression_dir + 'scran_normalization_regress_batch_' + regress_out_batch + '_individual_clustering_leiden_resolution_' + str(cluster_resolution) + '_cell_type_summary.txt'
# print_pseudobulk_clustering_mapping_cell_type_summary(adata, clustering_ct_summary_file, 'individual_leiden_clusters_' + str(cluster_resolution))



#######################
# Create list of ordered pseudobulk samples
# Ensure these samples pass some basic filters too
#######################
# Vector of length number of cells, where each element represents which cluster (ie pseudobulk sample) that cell gets assigned to
cluster_assignments = adata.obs['individual_leiden_clusters_' + str(cluster_resolution)]
# Extract orderd list of pseudobulk samples
min_cells_per_indi = 2500
ordered_pseudobulk_samples = extract_ordered_list_of_pseudobulk_samples(adata, genotyped_individuals_file, cluster_assignments, min_cells_per_indi)

#####################
# Get non-standardized expression
#####################
adata2 = adata.raw.to_adata()

'''
### NO LONGER USED!!!
sc.pp.highly_variable_genes(adata2)
# Get genes expressed in at least XX% of cells
min_fraction_of_cells = .01
sc.pp.filter_genes(adata2, min_cells=(adata2.X.shape[0])*min_fraction_of_cells)
# And genes that are highly variable
adata2 = adata2[:, adata2.var.highly_variable]
'''

# Make sure filtered genes line up (I know, this is all a little hacky..)
raw_ordered_genes = np.vstack((adata.var.index, adata.var[adata.var.columns[0]])).T
raw_ordered_genes2 = np.vstack((adata2.var.index, adata2.var[adata2.var.columns[0]])).T
if np.array_equal(raw_ordered_genes, raw_ordered_genes2) == False:
	print('assumption error')
	pdb.set_trace()

'''
#####################
# Get raw pseudobulk expression
#####################
if regress_out_batch == False:
	raw_pseudobulk_expression = generate_mean_cluster_pseudobulk_expression(adata2.X.toarray(), ordered_pseudobulk_samples, cluster_assignments, raw_ordered_genes)
else:
	raw_pseudobulk_expression = generate_mean_cluster_pseudobulk_expression(adata2.X, ordered_pseudobulk_samples, cluster_assignments, raw_ordered_genes)


#####################
# Save data to output
#####################	
gene_names_file = processed_pseudobulk_expression_dir + 'pseudobulk_scran_normalization_regress_batch_' + regress_out_batch + '_individual_clustering_leiden_resolution_' + str(cluster_resolution) + '_gene_names.txt'
np.savetxt(gene_names_file, raw_ordered_genes, fmt="%s", delimiter='\t')

# Sample names
sample_names_file = processed_pseudobulk_expression_dir + 'pseudobulk_scran_normalization_regress_batch_' + regress_out_batch + '_individual_clustering_leiden_resolution_' + str(cluster_resolution) + '_sample_names.txt'
np.savetxt(sample_names_file, ordered_pseudobulk_samples, fmt="%s", delimiter='\t')
# Generate pseudobulk covaraite file
'''
pseudobulk_covariate_file = processed_pseudobulk_expression_dir + 'pseudobulk_scran_normalization_regress_batch_' + regress_out_batch + '_individual_clustering_leiden_resolution_' + str(cluster_resolution) + '_sample_covariates.txt'
print_pseudobulk_covariate_file_from_cell_covariates(ordered_pseudobulk_samples, adata.obs, cluster_assignments, pseudobulk_covariate_file, donor_to_isg_score)

'''


#####################
# Normalize expression and generate expression pcs
#####################	
# Options for sample level normalization are currently 'none'
sample_level_normalization = 'qn'
# Options for gene level normalization are 'zscore' and 'ign'
gene_level_normalization = 'zscore'
# number of pcs
num_pcs = 200
# output root
pb_expression_output_root = processed_pseudobulk_expression_dir + 'pseudobulk_scran_normalization_regress_batch_' + regress_out_batch + '_individual_clustering_leiden_resolution_' + str(cluster_resolution) + '_' + sample_level_normalization + '_sample_norm_' + gene_level_normalization + '_gene_norm_'
normalize_expression_and_generate_expression_pcs(raw_pseudobulk_expression, sample_level_normalization, gene_level_normalization, num_pcs, pb_expression_output_root)


#####################
# Normalize expression and generate expression pcs
#####################	
# Options for sample level normalization are currently 'none'
sample_level_normalization = 'qn'
# Options for gene level normalization are 'zscore' and 'ign'
gene_level_normalization = 'ign'
# number of pcs
num_pcs = 200
# output root
pb_expression_output_root = processed_pseudobulk_expression_dir + 'pseudobulk_scran_normalization_regress_batch_' + regress_out_batch + '_individual_clustering_leiden_resolution_' + str(cluster_resolution) + '_' + sample_level_normalization + '_sample_norm_' + gene_level_normalization + '_gene_norm_'
normalize_expression_and_generate_expression_pcs(raw_pseudobulk_expression, sample_level_normalization, gene_level_normalization, num_pcs, pb_expression_output_root)


#####################
# Normalize expression and generate expression pcs
#####################	
# Options for sample level normalization are currently 'none'
sample_level_normalization = 'none'
# Options for gene level normalization are 'zscore' and 'ign'
gene_level_normalization = 'zscore'
# number of pcs
num_pcs = 200
# output root
pb_expression_output_root = processed_pseudobulk_expression_dir + 'pseudobulk_scran_normalization_regress_batch_' + regress_out_batch + '_individual_clustering_leiden_resolution_' + str(cluster_resolution) + '_' + sample_level_normalization + '_sample_norm_' + gene_level_normalization + '_gene_norm_'
normalize_expression_and_generate_expression_pcs(raw_pseudobulk_expression, sample_level_normalization, gene_level_normalization, num_pcs, pb_expression_output_root)


#####################
# Normalize expression and generate expression pcs
#####################	
# Options for sample level normalization are currently 'none'
sample_level_normalization = 'none'
# Options for gene level normalization are 'zscore' and 'ign'
gene_level_normalization = 'ign'
# number of pcs
num_pcs = 200
# output root
pb_expression_output_root = processed_pseudobulk_expression_dir + 'pseudobulk_scran_normalization_regress_batch_' + regress_out_batch + '_individual_clustering_leiden_resolution_' + str(cluster_resolution) + '_' + sample_level_normalization + '_sample_norm_' + gene_level_normalization + '_gene_norm_'
normalize_expression_and_generate_expression_pcs(raw_pseudobulk_expression, sample_level_normalization, gene_level_normalization, num_pcs, pb_expression_output_root)
'''
