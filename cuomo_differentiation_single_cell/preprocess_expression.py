import numpy as np 
import os
import sys
import pdb
import scanpy



def filter_expression_file_to_only_have_gene_symbols(normalized_expression_file, normalized_expression_gene_symbols_file):
	f = open(normalized_expression_file)
	t = open(normalized_expression_gene_symbols_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			num_header_features = len(data)
			t.write(line + '\n')
			continue
		# Error checking
		if len(data) != num_header_features:
			print('assumption eorororor')
			pdb.set_trace()
		new_gene = data[0].split('_')[1]
		t.write(new_gene + '\t' + '\t'.join(data[1:]) + '\n')
	t.close()
	f.close()


def get_dictionary_list_of_individuals_that_we_have_genotype_for(genotype_file):
	f = open(genotype_file)
	individuals = {}
	head_count = 0
	for line in f:
		if head_count == 0:
			line = line.rstrip()
			data = line.split()
			indis = data[1:]
			for indi in indis:
				individuals[indi] = 0
			head_count = head_count + 1
			continue
	f.close()
	return individuals

def make_cell_to_individual_mapping(meta_data_file, cell_to_individual_mapping_file, valid_individuals):
	f = open(meta_data_file)
	t = open(cell_to_individual_mapping_file, 'w')
	head_count = 0
	t.write('cell_name\tindividual_id\n')
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[85] in valid_individuals:
			valid_individuals[data[85]] = 1
			t.write(data[0] + '\t' + data[85] + '\n')
	f.close()
	t.close()

def recreate_cell_covariates(meta_data_file, recreated_cell_covariates_file, valid_individuals, genotype_pc_file, cell_state_file):
	# Extract mapping from indi id to first n genotype pcs
	n_pcs = 3
	mapping = {}
	f = open(genotype_pc_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 8:
			print('assumption error')
			pdb.set_trace()
		donor_id = data[0].split('"')[1]
		mapping[donor_id] = data[1:(1+n_pcs)]
	f.close()
	cell_state_mapping = {}
	f = open(cell_state_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			cell_state_names = np.asarray(data[1:])
			continue
		cell_id = data[0]
		cell_states = np.asarray(data[1:])
		# quick error checking
		if len(cell_states) != len(cell_state_names):
			print('assumption eororor')
		if cell_id in cell_state_mapping:
			print('assumption eroror')
			pdb.set_trace()
		cell_state_mapping[cell_id] = cell_states
	f.close()
	# Stream covariate file
	f = open(meta_data_file)
	t = open(recreated_cell_covariates_file, 'w')
	valid_cells = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write('cell_name_header\t' + '\t'.join(data) + '\tgenotype_pc1\tgenotype_pc2\tgenotype_pc3' + '\t' + '\t'.join(cell_state_names) + '\n')
			continue
		if data[85] in valid_individuals:
			genotype_pcs = mapping[data[85]]
			cell_states = cell_state_mapping[data[0]]
			t.write('\t'.join(data) + '\t' + '\t'.join(genotype_pcs) + '\t' + '\t'.join(cell_states) + '\n')
			valid_cells.append(True)
		else:
			valid_cells.append(False)
	f.close()
	t.close()
	return np.asarray(valid_cells)

def recreate_expression(normalized_expression_file, recreated_normalized_gene_expression_file, valid_cell_indices):
	f = open(normalized_expression_file)
	t = open(recreated_normalized_gene_expression_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = np.asarray(line.split())
		if head_count == 0:
			head_count = head_count + 1
			t.write('Gene_id\t' + '\t'.join(data[valid_cell_indices]) + '\n')
			continue
		gene_id = data[0]
		ensamble_id = gene_id.split('_')[0]
		counts = data[1:]
		t.write(gene_id + '\t' + '\t'.join(counts[valid_cell_indices]) + '\n')
	t.close()
	f.close()

def get_experiment_variable_from_cell_covariate_file(cell_covariates_file):
	f = open(cell_covariates_file)
	experiments = []
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			num_header_features = len(data)
			indices = np.where(np.asarray(data) == 'experiment')[0]
			if len(indices) != 1:
				print('assumption error')
				pdb.set_trace()
			index = indices[0]
			continue
		if len(data) != num_header_features:
			print('assumtpion eroror')
		experiments.append(data[index])
	return np.asarray(experiments)

def get_highly_variable_genes(recreated_normalized_gene_expression_gene_symbols_file, cell_covariates_file, method):
	experiment = get_experiment_variable_from_cell_covariate_file(cell_covariates_file)
	if method == 'scanpy_approach':
		adata = scanpy.read(recreated_normalized_gene_expression_gene_symbols_file).transpose()
		scanpy.pp.highly_variable_genes(adata)
		indices = adata.var['highly_variable']
	elif method == 'scran_approach':
		temp_expr = np.loadtxt(recreated_normalized_gene_expression_gene_symbols_file, dtype=str,delimiter='\t', comments='**')
		np.savetxt('experiment_temp.txt', experiment, fmt="%s", delimiter='\n')
		np.savetxt('expr_temp.txt', temp_expr[1:,1:], fmt="%s", delimiter='\t')
		os.system('Rscript get_scran_hvg.R expr_temp.txt experiment_temp.txt')
		num_genes = 1000
		scran_hvg_data = np.loadtxt('hvg.txt',dtype=str,delimiter='\t')
		hvg_bio = scran_hvg_data[1:,2].astype(float)
		indices = hvg_bio > sorted(hvg_bio,reverse=True)[num_genes]
	else:
		print('error: method ' + method + ' currently not implemented')
		pdb.set_trace()
	return indices

def filter_to_highly_variable_genes(recreated_normalized_gene_expression_file, recreated_normalized_hvg_gene_expression_file, highly_variable_gene_indices):
	f = open(recreated_normalized_gene_expression_file)
	t = open(recreated_normalized_hvg_gene_expression_file, 'w')
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count = head_count + 1
			t.write(line + '\n')
			num_header_features = len(data)
			continue
		# Error checks
		if num_header_features != len(data):
			print('assumptoin eroror')
			pdb.set_trace()
		if highly_variable_gene_indices[counter] == True:
			t.write(line + '\n')
		counter = counter + 1
	f.close()
	t.close()

def standardize_expression(normalized_expression_file, standardized_gene_expression_file):
	f = open(normalized_expression_file)
	t = open(standardized_gene_expression_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write('\t'.join(data) + '\n')
			continue
		counts = np.asarray(data[1:]).astype(float)
		standardized_counts = (counts - np.mean(counts))/np.std(counts)
		t.write(data[0] + '\t' + '\t'.join(standardized_counts.astype(str)) + '\n')
	f.close()
	t.close()

def standardize_expression_capped(normalized_expression_file, standardized_gene_expression_file, cap):
	f = open(normalized_expression_file)
	t = open(standardized_gene_expression_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write('\t'.join(data) + '\n')
			continue
		counts = np.asarray(data[1:]).astype(float)
		standardized_counts = (counts - np.mean(counts))/np.std(counts)
		standardized_counts[standardized_counts > cap] = cap
		standardized_counts[standardized_counts < -cap] = -cap
		t.write(data[0] + '\t' + '\t'.join(standardized_counts.astype(str)) + '\n')
	f.close()
	t.close()

# Generate expression PC loadings and variance explained of those expression PCs
def generate_pca_scores_and_variance_explained(filtered_standardized_sc_expression_file, num_pcs, filtered_cells_pca_file, filtered_cells_pca_ve_file):
	# Load in data
	X_full = np.loadtxt(filtered_standardized_sc_expression_file, dtype=str, delimiter='\t', comments='*')
	print(X_full[1:,1:].shape)
	X = np.transpose(X_full[1:,1:].astype(float))

	# Run PCA (via SVD)
	uuu, sss, vh = np.linalg.svd(np.transpose(X), full_matrices=False)
	svd_loadings = np.transpose(vh)[:,:num_pcs]

	# Save to output file
	np.savetxt(filtered_cells_pca_file, svd_loadings, fmt="%s", delimiter='\t')

	# Compute variance explained
	ve = (np.square(sss)/np.sum(np.square(sss)))[:num_pcs]
	np.savetxt(filtered_cells_pca_ve_file, ve, fmt="%s", delimiter='\n')



normalized_expression_file = sys.argv[1]
meta_data_file = sys.argv[2]
processed_genotype_dir = sys.argv[3]
gene_annotation_file = sys.argv[4]
genotype_pc_file = sys.argv[5]
cell_state_file = sys.argv[6]
processed_expression_dir = sys.argv[7]

genotype_file = processed_genotype_dir + 'genotype_missing_removed.txt'


##############################
# Reformat genotype data and get list of individuals that we have genotype data for
valid_individuals = get_dictionary_list_of_individuals_that_we_have_genotype_for(genotype_file)


###############################
# Create mapping from cell-id to individual id
# And filter cells to those for which we have genotype data for
cell_to_individual_mapping_file = processed_expression_dir + 'cell_individual_mapping.txt'
make_cell_to_individual_mapping(meta_data_file, cell_to_individual_mapping_file, valid_individuals)


###############################
# Re-create cell covariates
# And filter cells to those for which we have genotype data for
recreated_cell_covariates_file = processed_expression_dir + 'cell_covariates.txt'
valid_cell_indices = recreate_cell_covariates(meta_data_file, recreated_cell_covariates_file, valid_individuals, genotype_pc_file, cell_state_file)

###############################
# Re-create expression data
recreated_normalized_gene_expression_file = processed_expression_dir + 'normalized_expression_all_genotyped_cells.txt'
recreate_expression(normalized_expression_file, recreated_normalized_gene_expression_file, valid_cell_indices)


###############################
# Re-create expression data
recreated_normalized_gene_expression_gene_symbols_file = processed_expression_dir + 'normalized_expression_gene_symbols_all_genotyped_cells.txt'
filter_expression_file_to_only_have_gene_symbols(recreated_normalized_gene_expression_file, recreated_normalized_gene_expression_gene_symbols_file)


###############################
# Limit to highly variable genes
hvg_methods = ['scran', 'scanpy']

for hvg_method in hvg_methods:
	highly_variable_gene_indices = get_highly_variable_genes(recreated_normalized_gene_expression_gene_symbols_file, recreated_cell_covariates_file, hvg_method + '_approach')
	recreated_normalized_hvg_gene_expression_file = processed_expression_dir + 'normalized_expression_' + hvg_method + '_hvg_all_genotyped_cells.txt'
	filter_to_highly_variable_genes(recreated_normalized_gene_expression_file, recreated_normalized_hvg_gene_expression_file, highly_variable_gene_indices)


	###############################
	# Standardize expression data
	standardized_gene_expression_file = processed_expression_dir + 'standardized_normalized_expression_' + hvg_method + '_hvg_all_genotyped_cells.txt'
	standardize_expression(recreated_normalized_hvg_gene_expression_file, standardized_gene_expression_file)

	standardized_10_cap_gene_expression_file = processed_expression_dir + 'standardized_10_cap_normalized_expression_' + hvg_method + '_hvg_all_genotyped_cells.txt'
	standardize_expression_capped(recreated_normalized_hvg_gene_expression_file, standardized_10_cap_gene_expression_file, 10)

	###############################
	# Run PCA on standardized expression data
	num_pcs=200
	pca_loading_file = processed_expression_dir + 'standardized_normalized_expression_' + hvg_method + '_hvg_all_genotyped_cells_pca_loadings.txt'
	pca_pve_file = processed_expression_dir + 'standardized_normalized_expression_' + hvg_method + '_hvg_all_genotyped_cells_pca_pve.txt'
	generate_pca_scores_and_variance_explained(standardized_gene_expression_file, num_pcs, pca_loading_file, pca_pve_file)



	###############################
	# Run PCA on standardized expression data
	num_pcs=200
	pca_loading_file = processed_expression_dir + 'standardized_10_cap_normalized_expression_' + hvg_method + '_hvg_all_genotyped_cells_pca_loadings.txt'
	pca_pve_file = processed_expression_dir + 'standardized_10_cap_normalized_expression_' + hvg_method + '_hvg_all_genotyped_cells_pca_pve.txt'
	generate_pca_scores_and_variance_explained(standardized_10_cap_gene_expression_file, num_pcs, pca_loading_file, pca_pve_file)
