import numpy as np 
import os
import sys
import pdb







# Generate ordered list of cell types
def generate_ordered_list_of_cell_types(pseudobulk_covariate_file, cell_type_file):
	f = open(pseudobulk_covariate_file)
	cell_types = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		line_cell_type = '_'.join(data[1].split(' '))
		cell_types[line_cell_type] = 1
	f.close()
	t = open(cell_type_file, 'w')
	cell_type_arr = np.asarray(cell_types.keys())
	for cell_type in cell_type_arr:
		t.write(cell_type + '\n')
	t.close()
	return cell_type_arr

def filter_pseudobulk_covariate_file_to_single_cell_type(pseudobulk_covariate_file, cell_type, cell_type_sample_covariate_file):
	f = open(pseudobulk_covariate_file)
	t = open(cell_type_sample_covariate_file, 'w')
	head_count = 0
	indices = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		line_cell_type = '_'.join(data[1].split(' '))
		# Ignore samples not corresponding to this cell type
		if line_cell_type != cell_type:
			indices.append(False)
			continue
		indices.append(True)
		t.write(line + '\n')
	f.close()
	t.close()
	return np.asarray(indices)

# Standardize summed pseudobulk counts (samples X genes)
def standardize_pseudobulk_counts(raw_pseudobulk_expression, standardized_pseudobulk_expression_file):
	# Load in raw pseodobulk expression data
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

# For each cell type generate pseudobulk expression data
def generate_cell_type_pseudobulk_expression_data(cell_type, raw_pseudobulk_expression, pseudobulk_covariate_file, cell_type_raw_pseudobulk_expression_file, cell_type_pseudobulk_expression_file, cell_type_sample_covariate_file, cell_type_pca_loading_file, cell_type_pca_ve_file):
	# First filter pseudobulk covariate file to only samples from this cell type and return indices corresponding to this cell type
	cell_type_indices = filter_pseudobulk_covariate_file_to_single_cell_type(pseudobulk_covariate_file, cell_type, cell_type_sample_covariate_file)
	# Save cell type raw pseudobulk expression to output file
	cell_type_raw_pseudobulk_expression = raw_pseudobulk_expression[cell_type_indices, :]
	np.savetxt(cell_type_raw_pseudobulk_expression_file, cell_type_raw_pseudobulk_expression, fmt="%s", delimiter='\t')
	# Standardize summed pseudobulk counts (samples X genes)
	standardize_pseudobulk_counts(cell_type_raw_pseudobulk_expression, cell_type_pseudobulk_expression_file)
	# Compute PCs
	num_pcs=15
	generate_pca_scores_and_variance_explained(cell_type_pseudobulk_expression_file, num_pcs, cell_type_pca_loading_file, cell_type_pca_ve_file)

def pass_genotype_filter(genotype, thresh):
	booly = True
	rounded_geno = np.rint(genotype)
	num_samples = len(rounded_geno)

	num_zero = sum(rounded_geno==0)
	num_one = sum(rounded_geno == 1)
	num_two = sum(rounded_geno == 2)

	if num_zero + num_one + num_two != num_samples:
		print('assumption eroror')
		pdb.set_trace()

	fraction_zero = num_zero/float(num_samples)
	fraction_one = num_one/float(num_samples)
	fraction_two = num_two/float(num_samples)

	if fraction_zero < thresh:
		booly = False
	if fraction_one < thresh:
		booly = False
	if fraction_two < thresh:
		booly = False
	return booly

########################
# Step 1: Create file with all variant gene pairs such that gene is within $distanceKB of gene
########################
def extract_variant_gene_pairs_for_eqtl_testing(gene_file, gene_annotation_file, distance, genotype_data_dir, variant_gene_pair_file):
	# Extract gene list
	genes = np.loadtxt(gene_file, delimiter='\t',dtype=str)[:,0]
	# Get mapping from genes to (chrom_num, position)
	gene_mapping = get_mapping_from_gene_to_chromosome_position(gene_annotation_file, genes)
	# Open file handle to output file containing variant-gene pair tests and print header
	t = open(variant_gene_pair_file, 'w')
	t.write('Gene_id\tvariant_id\tchrom_num\tgene_tss\tvariant_position\n')

	# Fill in file containing lists of variant gene pairs for each chromosome iteratively
	used_genes = {}
	for chrom_num in range(1,23):
		print(chrom_num)
		# Create array where each element is a BP in this chromosome
		# 'Null' if no genes in distance BP of gene
		# Otherwise is a list of gene names
		chromosome = create_gene_chromsome(chrom_num, gene_mapping, distance)
		# Now loop through variants on this chromosome
		genotype_file = genotype_data_dir + 'clues_immvar_chrom_' + str(chrom_num) + '.DS.FORMAT'
		f = open(genotype_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split()
			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue
			if len(data) != 239:
				print('assumption error!')
			variant_id = data[0] + ':' + data[1]
			variant_chrom = int(data[0])
			variant_pos = int(data[1])
			# Simple error check
			if variant_chrom != chrom_num:
				print('assumption error')
				pdb.set_trace()
			# No genes within 10KB of variant
			if chromosome[variant_pos] == 'Null':
				continue
			genotype = np.asarray(data[2:]).astype(float)
			maf = get_maf(genotype)
			if maf < .1:
				print('skipped variant' + '\t' + str(maf))
				continue
			if pass_genotype_filter(genotype, .05) == False:
				continue
			# List of genes that variant maps to
			mapping_genes = chromosome[variant_pos].split(':')
			if len(mapping_genes) != len(np.unique(mapping_genes)):
				print('assumption erroro')
				pdb.set_trace()
			for gene_id in mapping_genes:
				# THIS IS A VARIANT-GENE PAIR WE WILL TEST
				# PRINT TO OUTPUT
				used_genes[gene_id] = 1
				t.write(gene_id + '\t' + variant_id + '\t' + str(chrom_num) + '\t' + str(gene_mapping[gene_id][1]) + '\t' + str(variant_pos) + '\n')
		f.close()
	t.close()

# Get mapping from genes to (chrom_num, position)
def get_mapping_from_gene_to_chromosome_position(gene_annotation_file, genes):
	# Convert gene array to dictionary
	gene_mapping = {}
	for gene in genes:
		gene_mapping[gene] = [0,0]
	# Fill in gene positions
	f = open(gene_annotation_file)
	used_genes = {}
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
		if gene_id not in gene_mapping:
			continue
		used_genes[gene_id] = 1
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
		# Extract relevent info on gene: Chrom num and TSS
		chrom_num = int(data[1].split('hr')[1])
		strand = data[4]
		gene_id = data[5]
		if strand == '+':
			tss = start
		elif strand == '-':
			tss = end
		else:
			print('assumption error while processing gene annotation file')
		# Add info to gene dictionary
		if gene_mapping[gene_id][0] != 0:
			if gene_mapping[gene_id][0] != chrom_num and gene_mapping[gene_id][1] != tss:
				gene_mapping.pop(gene_id)
		else:
			gene_mapping[gene_id] = (chrom_num, tss)
	f.close()
	return gene_mapping

# Create array where each element is a BP in this chromosome
# 'Null' if no genes in distance BP of gene
# Otherwise is a list of gene names
def create_gene_chromsome(chrom_num, gene_mapping, distance):
	# Initialize chromosome
	chromosome = ['Null']*250000000
	# Loop through genes
	for gene_id in gene_mapping.keys():
		# Extract dictionary on position of gene
		gene_tuple = gene_mapping[gene_id]
		gene_chrom_num = gene_tuple[0]
		gene_tss = gene_tuple[1]
		# Make sure gene is on correct chromosome
		if gene_chrom_num != chrom_num:
			continue
		# Iterate across positions in a distance KB window around gene and add gene to chromosome
		for pos in range((gene_tss-distance),(gene_tss+distance+1)):
			if chromosome[pos] == 'Null':
				chromosome[pos] = gene_id
			else:
				chromosome[pos] = chromosome[pos] + ':' + gene_id
	return chromosome

def get_maf(genotype):
	af = np.sum(genotype)/(2.0*len(genotype))
	if af > .5:
		maf = 1.0 - af
	else:
		maf = af
	if maf > .5 or maf < 0.0:
		print('assumption error')
		pdb.set_trace()
	return maf
def get_ordered_list_of_gene_names(ld_pruned_variant_gene_pair_file):
	gene_names = []
	f = open(ld_pruned_variant_gene_pair_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_names.append(data[0])
	return np.asarray(gene_names)

def get_ordered_list_of_variant_names(ld_pruned_variant_gene_pair_file):
	variant_names = []
	dicti = {}
	f = open(ld_pruned_variant_gene_pair_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_names.append(data[1])
		dicti[data[1]] = 1
	return np.asarray(variant_names), dicti

def create_mapping_from_gene_names_to_expression_vectors(sc_expression_file, gene_names_file):
	# Get gene names
	gene_names = np.loadtxt(gene_names_file,dtype=str, delimiter='\t')[:,0]
	# Load in expression matrix (every column corresponds to a gene)
	expression_matrix = np.loadtxt(sc_expression_file, dtype=str, delimiter='\t')
	#np.save('temp.npy', expression_matrix)
	#expression_matrix = np.load('temp.npy')
	#print('loaded')
	# Create mapping
	mapping = {}
	for index, gene_name in enumerate(gene_names):
		mapping[gene_name] = expression_matrix[:, index]
	del expression_matrix
	return mapping

def create_mapping_from_gene_names_to_expression_vector_float(sc_expression_file, gene_names_file):
	# Get gene names
	gene_names = np.loadtxt(gene_names_file,dtype=str, delimiter='\t')[:,0]
	# Load in expression matrix (every column corresponds to a gene)
	expression_matrix = np.loadtxt(sc_expression_file, dtype=str, delimiter='\t')
	# Create mapping
	mapping = {}
	for index, gene_name in enumerate(gene_names):
		mapping[gene_name] = np.asarray(expression_matrix[:, index]).astype(float)
	return mapping

def generate_single_cell_expression_eqtl_training_data(ld_pruned_variant_gene_pair_file, sc_expression_file, gene_names_file, single_cell_expression_eqtl_traing_data_file):
	# Get ordered list of gene names (this will be the order that the output file will be saved in)
	ordered_gene_names = get_ordered_list_of_gene_names(ld_pruned_variant_gene_pair_file)
	# Create mapping from gene names to expression vectors
	gene_name_to_expression_vector = create_mapping_from_gene_names_to_expression_vectors(sc_expression_file, gene_names_file)
	# Now convert vectors to strings
	for gene_name in gene_name_to_expression_vector.keys():
		gene_name_to_expression_vector[gene_name] = '\t'.join(gene_name_to_expression_vector[gene_name])
	# print to output file
	t = open(single_cell_expression_eqtl_traing_data_file, 'w')
	for gene_name in ordered_gene_names:
		# Use map to get expression vector corresponding to this gene
		if gene_name not in gene_name_to_expression_vector:
			pdb.set_trace()
		expression_string = gene_name_to_expression_vector[gene_name]
		# Print to output file
		t.write(expression_string + '\n')
	t.close()


def get_cell_level_ordered_individaul_array(cell_level_info_file):
	array = []
	head_count = 0
	f = open(cell_level_info_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		array.append(data[1])
	f.close()
	return np.asarray(array)


def get_genotype_level_ordered_individual_array(genotype_data_dir):
	chrom_num = 1
	# genotype_file = genotype_data_dir + 'chr' + str(chrom_num) + '.genotypes.matrix.eqtl.txt'
	genotype_file = genotype_data_dir + 'clues_immvar_chrom_' + str(chrom_num) + '.DS.FORMAT'
	f = open(genotype_file)
	head_count = 0
	for line in f:
		if head_count == 0:
			line = line.rstrip()
			data = line.split()
			indi = data[2:]
			head_count = head_count + 1
			continue
	f.close()
	return np.asarray(indi)

# Create mapping from variants in variat_list to genotype vectors
def create_mapping_from_variants_to_genotype(variant_list, genotype_data_dir):
	variants = {}
	for chrom_num in range(1,23):
		# genotype_file = genotype_data_dir + 'chr' + str(chrom_num) + '.genotypes.matrix.eqtl.txt'
		genotype_file = genotype_data_dir + 'clues_immvar_chrom_' + str(chrom_num) + '.DS.FORMAT'
		f = open(genotype_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split()
			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue
			# Simple error checking
			if len(data) != 239:
				print('assumption error!')
			variant_id = data[0] + ':' + data[1]
			# Limit to variants in variant_list
			if variant_id not in variant_list:
				continue
			variants[variant_id] = np.asarray(data[2:])
		f.close()
	if len(variants) != len(variant_list):
		print('assumption error')
		pdb.set_trace()
	return variants

def construct_genotype_matrix(ld_pruned_variant_gene_pair_file, genotype_data_dir, cell_level_info_file, single_cell_genotype_eqtl_training_data_file):
	variant_names, variant_list = get_ordered_list_of_variant_names(ld_pruned_variant_gene_pair_file)

	ordered_individuals_cell_level = get_cell_level_ordered_individaul_array(cell_level_info_file)
	
	ordered_individuals_genotype_level = get_genotype_level_ordered_individual_array(genotype_data_dir)

	# Create mapping from variant_id to genotype vectors
	variant_to_genotype = create_mapping_from_variants_to_genotype(variant_list, genotype_data_dir)

	# Create array mapping from genotype-level array to single cell-level array
	mapping_array = []
	converter = {}
	for i, indi in enumerate(ordered_individuals_genotype_level):
		converter[indi] = i 
	for indi in ordered_individuals_cell_level:
		mapping_array.append(converter[indi])
	mapping_array = np.asarray(mapping_array)
	# Simple error checking
	if np.array_equal(ordered_individuals_cell_level, ordered_individuals_genotype_level[mapping_array]) == False:
		print('assumption error')
		pdb.set_trace()

	# print new genotype file
	t = open(single_cell_genotype_eqtl_training_data_file, 'w')
	for variant_id in variant_names:
		genotype_vector = (variant_to_genotype[variant_id])
		cell_level_genotype_vector = genotype_vector[mapping_array]
		t.write('\t'.join(cell_level_genotype_vector) + '\n')
	t.close()

def construct_sample_overlap_file(cell_type_sc_sample_covariate_file, cell_type_sample_overlap_file):
	ordered_individuals_cell_level = get_cell_level_ordered_individaul_array(cell_type_sc_sample_covariate_file)
	unique_indis = np.unique(ordered_individuals_cell_level)
	mapping_from_cell_to_number = {}
	for i,indi in enumerate(unique_indis):
		mapping_from_cell_to_number[indi] = i
	t = open(cell_type_sample_overlap_file, 'w')
	for indi in ordered_individuals_cell_level:
		t.write(str(mapping_from_cell_to_number[indi]) + '\n')
	t.close() 


def generate_test_level_covariate_file(cell_type_eqtl_variant_gene_pairs_file, cell_type_sc_raw_expression_file, gene_id_file, cell_type_test_info_file):
	# Get ordered list of gene names (this will be the order that the output file will be saved in)
	ordered_gene_names = get_ordered_list_of_gene_names(cell_type_eqtl_variant_gene_pairs_file)
	# Create mapping from gene names to expression vectors
	gene_name_to_expression_vector = create_mapping_from_gene_names_to_expression_vector_float(cell_type_sc_raw_expression_file, gene_id_file)
	# print to output file
	t = open(cell_type_test_info_file, 'w')
	t.write('total_counts\tpercent_expressed_cells\n')
	for gene_name in ordered_gene_names:
		# Use map to get expression vector corresponding to this gene
		expression_vector = gene_name_to_expression_vector[gene_name]
		total_counts = np.sum(expression_vector)
		fraction_expressed = float(np.sum(expression_vector != 0.0))/len(expression_vector)
		# Print to output file
		t.write(str(total_counts) + '\t' + str(fraction_expressed) + '\n')
	t.close()

def construct_library_size_file(cell_type_sc_sample_covariate_file, cell_type_library_size_file):
	f = open(cell_type_sc_sample_covariate_file)
	t = open(cell_type_library_size_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		t.write(data[9] + '\n')
	f.close()
	t.close()

def construct_covariate_file(sc_sample_covariate_file, eqtl_covariate_file):
	# First need to extract genotype pcs
	f = open(sc_sample_covariate_file)
	head_count = 0
	geno_pcs = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			# Get indices that contain genotype pcs
			geno_pc_indices = []
			for index, column_name in enumerate(data):
				if column_name.startswith('genotype_pc'):
					geno_pc_indices.append(index)
			geno_pc_indices = np.asarray(geno_pc_indices)
			continue
		line_geno_pcs = np.asarray(data)[geno_pc_indices].astype(float)
		geno_pcs.append(line_geno_pcs)
	f.close()
	geno_pcs = np.asarray(geno_pcs)
	# Merge geno pcs and expression pcs
	#covs = np.hstack((geno_pcs, pcs))
	np.savetxt(eqtl_covariate_file, geno_pcs, fmt="%s", delimiter='\t')

def generate_latent_factor_interaction_eqtl_input_files(genotype_data_dir, expression_file, sample_covariate_file, expression_pcs_file, eqtl_variant_gene_pairs_file, eqtl_expression_file, eqtl_genotype_file, distance, gene_annotation_file, gene_id_file, eqtl_sample_overlap_file, eqtl_covariate_file, eqtl_lf_file, num_pcs):
	########################
	# Step 1: Create file with all variant gene pairs such that gene is within $distanceKB of gene
	########################
	extract_variant_gene_pairs_for_eqtl_testing(gene_id_file, gene_annotation_file, distance, genotype_data_dir, eqtl_variant_gene_pairs_file)

	########################
	# Step 2: Generate expression matrix
	########################
	print('expr')
	generate_single_cell_expression_eqtl_training_data(eqtl_variant_gene_pairs_file, expression_file, gene_id_file, eqtl_expression_file)

	########################
	# Step 3: Generate Genotype matrix
	########################
	print('geno')
	construct_genotype_matrix(eqtl_variant_gene_pairs_file, genotype_data_dir, sample_covariate_file, eqtl_genotype_file)

	########################
	# Step 4: Generate sample overlap file
	########################
	print('sample overlap')
	construct_sample_overlap_file(sample_covariate_file, eqtl_sample_overlap_file)

	########################
	# Step 5: Generate covariate file
	########################
	construct_covariate_file(sample_covariate_file, eqtl_covariate_file)


	########################
	# Step 6: Generate latent factor file
	########################
	pcs_full = np.loadtxt(expression_pcs_file)
	pcs = pcs_full[:, :num_pcs]
	np.savetxt(eqtl_lf_file, pcs, fmt="%s", delimiter='\t')


###################
# Command Line args
###################
processed_expression_dir = sys.argv[1]  # Input dir
genotype_data_dir = sys.argv[2]  # Input dir
gene_annotation_file = sys.argv[3]
latent_factor_interaction_eqtl_dir = sys.argv[4]  # Output dir


###################
# Input files
###################
# Pseudobulk gene names 
gene_id_file = processed_expression_dir + 'cluster_tmm_ign_pseudobulk_leiden_no_cap_2.5_gene_names.txt'


###################
# Extract gene ids
###################
gene_id_mat = np.loadtxt(gene_id_file, dtype=str)
gene_symbol_ids = gene_id_mat[:,0]
ensamble_ids = gene_id_mat[:,1]
# Simple error checking
if len(np.unique(ensamble_ids)) != len(np.unique(gene_symbol_ids)):
	print('assumption error')
	pdb.set_trace()
if len(np.unique(ensamble_ids)) != len(ensamble_ids):
	print('assumption error')
	pdb.set_trace()



###################
# Extract latent factor eqtl input data
###################
# Max distance between variant and gene tss
distance=10000
num_pcs = 50

# Input files
expression_file = processed_expression_dir + 'cluster_tmm_ign_pseudobulk_leiden_no_cap_2.5_log_tmm_normalized_expression.txt'
sample_covariate_file = genotype_data_dir + 'pseudobulk_sample_covariates_with_genotype_pcs.txt'
expression_pcs_file = processed_expression_dir + 'cluster_tmm_ign_pseudobulk_leiden_no_cap_2.5_pca_scores.txt'

# Output files
eqtl_variant_gene_pairs_file = latent_factor_interaction_eqtl_dir + 'latent_factor_interaction_no_cap_2.5_eqtl_input_variant_gene_pairs.txt'
eqtl_expression_file = latent_factor_interaction_eqtl_dir + 'latent_factor_interaction_no_cap_2.5_eqtl_input_expression.txt'
eqtl_genotype_file = latent_factor_interaction_eqtl_dir + 'latent_factor_interaction_no_cap_2.5_eqtl_input_genotype.txt'
eqtl_sample_overlap_file = latent_factor_interaction_eqtl_dir + 'latent_factor_interaction_no_cap_2.5_eqtl_input_sample_overlap.txt'
eqtl_covariate_file = latent_factor_interaction_eqtl_dir + 'latent_factor_interaction_no_cap_2.5_eqtl_input_covariates.txt'
eqtl_lf_file = latent_factor_interaction_eqtl_dir + 'latent_factor_interaction_no_cap_2.5_eqtl_input_latent_factors.txt'
generate_latent_factor_interaction_eqtl_input_files(genotype_data_dir, expression_file, sample_covariate_file, expression_pcs_file, eqtl_variant_gene_pairs_file, eqtl_expression_file, eqtl_genotype_file, distance, gene_annotation_file, gene_id_file, eqtl_sample_overlap_file, eqtl_covariate_file, eqtl_lf_file, num_pcs)



