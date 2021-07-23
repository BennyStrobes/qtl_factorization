import numpy as np 
import os
import sys
import pdb

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
		ensamble_id = data[0].split('.')[0]
		if ensamble_id not in gene_mapping:
			continue
		chrom_num = data[1]
		if chrom_num == 'chrX' or chrom_num == 'chrY' or chrom_num == 'chrM':
			continue
		# Simple error checking
		start = int(data[2])
		end = int(data[3])
		if start > end:
			print('assumption eroror')
			pdb.set_trace()
		# Extract relevent info on gene: Chrom num and TSS
		chrom_num = int(data[1].split('hr')[1])
		strand = data[4]
		if strand == '+':
			tss = start
		elif strand == '-':
			tss = end
		else:
			print('assumption error while processing gene annotation file')
		# Add info to gene dictionary
		if gene_mapping[ensamble_id][0] != 0:
			print('miss')
			pdb.set_trace()
			if gene_mapping[ensamble_id][0] != chrom_num and gene_mapping[ensamble_id][1] != tss:
				gene_mapping.pop(ensamble_id)
		else:
			gene_mapping[ensamble_id] = (chrom_num, tss)
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



########################
# Step 1: Create file with all variant gene pairs such that gene is within $distanceKB of gene
########################
def extract_variant_gene_pairs_for_eqtl_testing(gene_file, gene_annotation_file, distance, genotype_file, variant_gene_pair_file):
	# Extract gene list
	genes_temp = np.loadtxt(gene_file, delimiter='\t',dtype=str, comments='**')[1:,0]
	genes = []
	for gene in genes_temp:
		genes.append(gene.split('_')[0])
	genes = np.asarray(genes)
	# Get mapping from genes to (chrom_num, position)
	gene_mapping = get_mapping_from_gene_to_chromosome_position(gene_annotation_file, genes)
	# Open file handle to output file containing variant-gene pair tests and print header
	t = open(variant_gene_pair_file, 'w')
	t.write('Gene_id\tvariant_id\tchrom_num\tgene_tss\tvariant_position\n')
	# Fill in file containing lists of variant gene pairs for each chromosome iteratively
	for chrom_num in range(1,23):
		print(chrom_num)
		# Create array where each element is a BP in this chromosome
		# 'Null' if no genes in distance BP of gene
		# Otherwise is a list of gene names
		chromosome = create_gene_chromsome(chrom_num, gene_mapping, distance)
		# Now loop through variants on this chromosome
		f = open(genotype_file)
		skipped = 0
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue
			if len(data) != 106:
				print('assumption error!')
				pdb.set_trace()
			variant_id = data[0]
			variant_chrom = int(variant_id.split(':')[0])
			variant_pos = int(variant_id.split(':')[1])
			# Ignore variants from other chromosomes
			if variant_chrom != chrom_num:
				continue
			# No genes within distance of variant
			if chromosome[variant_pos] == 'Null':
				continue
			# List of genes that variant maps to
			mapping_genes = chromosome[variant_pos].split(':')
			for gene_id in mapping_genes:
				# THIS IS A VARIANT-GENE PAIR WE WILL TEST
				# PRINT TO OUTPUT
				t.write(gene_id + '\t' + variant_id + '\t' + str(chrom_num) + '\t' + str(gene_mapping[gene_id][1]) + '\t' + str(variant_pos) + '\n')
		f.close()
	t.close()

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

def create_mapping_from_gene_names_to_expression_vectors(sc_expression_file):
	# Extract gene list
	genes_full = np.loadtxt(sc_expression_file, delimiter='\t',dtype=str, comments='**')
	genes_temp = genes_full[1:, 0]
	genes = []
	for gene in genes_temp:
		genes.append(gene.split('_')[0])
	gene_names = np.asarray(genes)
	expression_matrix = np.transpose(genes_full[1:,1:])

	# Create mapping
	mapping = {}
	for index, gene_name in enumerate(gene_names):
		mapping[gene_name] = expression_matrix[:, index]
	return mapping



def generate_single_cell_expression_eqtl_training_data(ld_pruned_variant_gene_pair_file, sc_expression_file, single_cell_expression_eqtl_traing_data_file):
	# Get ordered list of gene names (this will be the order that the output file will be saved in)
	ordered_gene_names = get_ordered_list_of_gene_names(ld_pruned_variant_gene_pair_file)
	# Create mapping from gene names to expression vectors
	gene_name_to_expression_vector = create_mapping_from_gene_names_to_expression_vectors(sc_expression_file)
	# print to output file
	t = open(single_cell_expression_eqtl_traing_data_file, 'w')
	for gene_name in ordered_gene_names:
		# Use map to get expression vector corresponding to this gene
		if gene_name not in gene_name_to_expression_vector:
			pdb.set_trace()
		expression_vector = gene_name_to_expression_vector[gene_name]
		# Print to output file
		t.write('\t'.join(expression_vector) + '\n')
	t.close()

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
		array.append(data[85])
	f.close()
	return np.asarray(array)

def get_genotype_level_ordered_individual_array(genotype_file):
	f = open(genotype_file)
	head_count = 0
	for line in f:
		if head_count == 0:
			line = line.rstrip()
			data = line.split()
			indi = data[1:]
			head_count = head_count + 1
			continue
	f.close()
	return np.asarray(indi)

# Create mapping from variants in variat_list to genotype vectors
def create_mapping_from_variants_to_genotype(variant_list, genotype_file):
	variants = {}
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
		if len(data) != 106:
			print('assumption error!')
		variant_id = data[0]
		# Limit to variants in variant_list
		if variant_id not in variant_list:
			continue
		variants[variant_id] = np.asarray(data[1:]).astype(float)
	f.close()
	if len(variants) != len(variant_list):
		print('assumption error')
		pdb.set_trace()
	return variants



def construct_genotype_matrix(ld_pruned_variant_gene_pair_file, genotype_file, cell_level_info_file, single_cell_genotype_eqtl_training_data_file):
	variant_names, variant_list = get_ordered_list_of_variant_names(ld_pruned_variant_gene_pair_file)

	ordered_individuals_cell_level = get_cell_level_ordered_individaul_array(cell_level_info_file)

	ordered_individuals_genotype_level = get_genotype_level_ordered_individual_array(genotype_file)

	# Create mapping from variant_id to genotype vectors
	variant_to_genotype = create_mapping_from_variants_to_genotype(variant_list, genotype_file)

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
		genotype_vector = (variant_to_genotype[variant_id]).astype(str)
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

def construct_covariate_file(cell_info_file, pc_file, num_pcs, eqtl_covariate_file):
	pcs = np.loadtxt(pc_file, dtype=str, delimiter='\t')[:, :num_pcs]
	geno_pcs = []
	head_count = 0
	f = open(cell_info_file)
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			num_headers = len(data)
			continue
		if len(data) != num_headers:
			print('assumption erooror')
			pdb.set_trace()
		geno_pcs.append(data[94:97])
	f.close()
	geno_pcs = np.asarray(geno_pcs)
	cov = np.hstack((pcs, geno_pcs))
	np.savetxt(eqtl_covariate_file, cov, fmt="%s", delimiter='\t')

def generate_eqtl_input_files(distance, gene_annotation_file, genotype_file, expression_file, cell_info_file, pc_file, eqtl_variant_gene_pairs_file, eqtl_expression_file, eqtl_genotype_file, eqtl_sample_overlap_file, eqtl_covariate_file):
	########################
	# Step 1: Create file with all variant gene pairs such that gene is within $distanceKB of gene
	########################
	extract_variant_gene_pairs_for_eqtl_testing(expression_file, gene_annotation_file, distance, genotype_file, eqtl_variant_gene_pairs_file)

	########################
	# Step 2: Generate expression matrix
	########################
	generate_single_cell_expression_eqtl_training_data(eqtl_variant_gene_pairs_file, expression_file, eqtl_expression_file)
	
	########################
	# Step 3: Generate Genotype matrix
	########################
	construct_genotype_matrix(eqtl_variant_gene_pairs_file, genotype_file, cell_info_file, eqtl_genotype_file)

	########################
	# Step 4: Generate sample overlap file
	########################
	construct_sample_overlap_file(cell_info_file, eqtl_sample_overlap_file)

	########################
	# Step 5: Generate covariate file
	########################
	construct_covariate_file(cell_info_file, pc_file, 50, eqtl_covariate_file)


## Command line arguments
processed_expression_dir = sys.argv[1]
processed_genotype_dir = sys.argv[2]
gene_annotation_file = sys.argv[3]
standard_eqtl_input_dir = sys.argv[4]


# Input files
genotype_file = processed_genotype_dir + 'genotype_missing_removed.txt'
expression_file = processed_expression_dir + 'standardized_10_cap_normalized_expression_scran_hvg_all_genotyped_cells.txt'
cell_info_file = processed_expression_dir + 'cell_covariates.txt'
pc_file = processed_expression_dir + 'standardized_10_cap_normalized_expression_scran_hvg_all_genotyped_cells_pca_loadings.txt'

# Output files
output_prefix = 'scran_1000_hvg_'
eqtl_variant_gene_pairs_file = standard_eqtl_input_dir + output_prefix + 'eqtl_input_variant_gene_pairs.txt'
eqtl_expression_file = standard_eqtl_input_dir + output_prefix + 'eqtl_input_expression.txt'
eqtl_genotype_file = standard_eqtl_input_dir + output_prefix + 'eqtl_input_genotype.txt'
eqtl_sample_overlap_file = standard_eqtl_input_dir + output_prefix + 'eqtl_input_sample_overlap.txt'
eqtl_covariate_file = standard_eqtl_input_dir + output_prefix + 'eqtl_input_covariates.txt'

distance=100000

generate_eqtl_input_files(distance, gene_annotation_file, genotype_file, expression_file, cell_info_file, pc_file, eqtl_variant_gene_pairs_file, eqtl_expression_file, eqtl_genotype_file, eqtl_sample_overlap_file, eqtl_covariate_file)


