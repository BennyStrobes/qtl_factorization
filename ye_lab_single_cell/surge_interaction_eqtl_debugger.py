import numpy as np 
import os
import sys
import pdb


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
	genotype_file = genotype_data_dir + 'clues_immvar_chrom_' + str(chrom_num) + '.DS2.FORMAT'
	f = open(genotype_file)
	head_count = 0
	for line in f:
		if head_count == 0:
			line = line.rstrip()
			data = line.split()
			indi = data[3:]
			head_count = head_count + 1
			continue
	f.close()
	return np.asarray(indi)





# Create mapping from variants in variat_list to genotype vectors
def create_mapping_from_variants_to_genotype(variant_list, genotype_data_dir):
	variants = {}
	for chrom_num in range(1,23):
		# genotype_file = genotype_data_dir + 'chr' + str(chrom_num) + '.genotypes.matrix.eqtl.txt'
		genotype_file = genotype_data_dir + 'clues_immvar_chrom_' + str(chrom_num) + '.DS2.FORMAT'
		f = open(genotype_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split()
			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue
			variant_id = data[2]
			# Limit to variants in variant_list
			if variant_id not in variant_list:
				continue
			variants[variant_id] = np.asarray(data[3:])
		f.close()
	if len(variants) != len(variant_list):
		print('assumption error')
		pdb.set_trace()
	return variants



output_stem = sys.argv[1]


lf_num = 0

full_results_file = output_stem + 'interaction_only_interaction_eqtl_results_latent_factor_' + str(lf_num+1) + '_merged.txt'
latent_factor_file = output_stem + 'interaction_only_surge_latent_factors.txt'


lf_mat = np.loadtxt(latent_factor_file)
lf_vec = lf_mat[:, lf_num]


genes = {}
variant_gene_pairs = {}
variants = {}


f = open(full_results_file)
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	#lf_p_vec = np.asarray(data[3:]).astype(float)
	#lf_p = lf_p_vec[lf_num]
	variant_id = data[0]
	gene_id = data[1]
	lf_p = float(data[4])
	if lf_p < 1e-200:
		pdb.set_trace()
		genes[gene_id] = 1
		variants[variant_id] = 1
		variant_gene_pairs[variant_id + '_' + gene_id] = 1
f.close()

#variant_to_geno = {}
genotype_data_dir = '/scratch16/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/processed_genotype/'

variant_to_geno = create_mapping_from_variants_to_genotype(variants, genotype_data_dir)


cell_level_info_file = genotype_data_dir +'pseudobulk_sample_covariates_with_genotype_pcs.txt'

ordered_individuals_cell_level = get_cell_level_ordered_individaul_array(cell_level_info_file)
ordered_individuals_genotype_level = get_genotype_level_ordered_individual_array(genotype_data_dir)

# Create array mapping from genotype-level array to single cell-level array
mapping_array = []
converter = {}
for i, indi in enumerate(ordered_individuals_genotype_level):
	converter[indi] = i 
for indi in ordered_individuals_cell_level:
	mapping_array.append(converter[indi])
mapping_array = np.asarray(mapping_array)

arr = []
for variant_id in variants.keys():
	ind_genotype = variant_to_geno[variant_id].astype(float)
	samp_genotype = ind_genotype[mapping_array]
	stdev = np.std((samp_genotype-np.mean(samp_genotype)))
	arr.append(stdev)
pdb.set_trace()
