import numpy as np 
import os
import sys
import pdb



def extract_significant_genes(sig_eqtl_file, chrom_num):
	gene_mapping = {}
	f = open(sig_eqtl_file)
	head_count = 0
	chrom_string = str(chrom_num) + ':'
	for line in f:
		line = line.rstrip()
		if line.startswith(chrom_string) == False:
			continue
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Error checking
		# Extract relevent fields
		gene_name = data[1]
		sentinal_snp = data[0]
		# more error checking
		if gene_name in gene_mapping:
			print('assumption errorror')
			pdb.set_trace()
		gene_mapping[gene_name] = []
	f.close()
	return gene_mapping


def map_genes_to_all_variants(all_eqtl_file, chrom_num, gene_mapping):
	f = open(all_eqtl_file)
	chrom_string = str(chrom_num) + ':'
	for line in f:
		line = line.rstrip()
		if line.startswith(chrom_string) == False:
			continue
		data = line.split('\t')
		if len(data) != 5:
			print('assumption erorro')
			pdb.set_trace()
		gene_name = data[1]
		if gene_name not in gene_mapping:
			continue
		variant_id = data[0]
		variant_info = variant_id.split(':')
		variant_id_alt = variant_info[0] + ':' + variant_info[1] + ':' + variant_info[3] + ':' + variant_info[2]
		beta = data[2]
		std_err_beta = data[3]
		var_beta = str(np.square(float(std_err_beta)))
		pvalue = data[4]
		gene_mapping[gene_name].append((variant_id, beta, var_beta, pvalue))
	f.close()
	return gene_mapping

def make_variant_to_sample_genotype_vec_mapping(genotype_file):
	f = open(genotype_file)
	mapping = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			continue
		# Simple error checking
		if len(data) != 240:
			print('assumption error!')
		variant_id = data[2]
		genotype_vec = np.asarray(data[3:]).astype(float)
		# Simple error checking
		if variant_id in mapping:
			print('repeat asssumption errorr')
			pdb.set_trace()
		mapping[variant_id] = genotype_vec
	f.close()
	return mapping

def reorganize_data_for_susie_input(eqtl_variants, variant_to_genotype_vec_mapping, loading_vec, donor_genotype_to_cell_genotype_mapping):
	variant_names = []
	snp_positions = []
	eqtl_betas = []
	eqtl_beta_vars = []
	eqtl_pvalues = []
	eqtl_genotype_mat = []
	for eqtl_variant_tuple in eqtl_variants:
		eqtl_variant = eqtl_variant_tuple[0]
		eqtl_variant_beta = float(eqtl_variant_tuple[1])
		eqtl_variant_beta_var = float(eqtl_variant_tuple[2])
		eqtl_variant_pvalue = float(eqtl_variant_tuple[3])
		# Save data to output vecs
		variant_names.append(eqtl_variant)
		snp_positions.append(eqtl_variant.split(':')[1])
		eqtl_betas.append(eqtl_variant_beta)
		eqtl_beta_vars.append(eqtl_variant_beta_var)
		eqtl_pvalues.append(eqtl_variant_pvalue)
		# convert from donor genotype to cell level genotype
		donor_genotype = variant_to_genotype_vec_mapping[eqtl_variant]
		cell_genotype = donor_genotype[donor_genotype_to_cell_genotype_mapping]
		eqtl_genotype_mat.append(cell_genotype*loading_vec)

	# Put into organized matrix (susie format)
	eqtl_mat = np.transpose(np.vstack((np.asarray(variant_names), np.asarray(snp_positions), np.asarray(eqtl_betas).astype(str), np.asarray(eqtl_beta_vars).astype(str), np.asarray(eqtl_pvalues).astype(str))))
	# Add headers
	eqtl_mat_w_header = np.vstack((['snp', 'position', 'beta', 'varbeta', 'pvalue'], eqtl_mat))

	# Create matrix of genotype for this gene of dimension num_snpsXnum_individuals
	eqtl_genotype_mat = np.asarray(eqtl_genotype_mat)
	# Generate ld matrix of dimension num_snpsXnum_snps
	ld_mat = np.corrcoef(eqtl_genotype_mat)
	return True, eqtl_mat_w_header, ld_mat

def get_donor_genotype_to_cell_genotype_index_mapping(genotype_file, sample_names_file):
	f = open(genotype_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			donor_names = np.asarray(data[3:])
			donor_mapping = {}
			for i, donor_name in enumerate(donor_names):
				donor_mapping[donor_name] = i
			continue
		break
	f.close()
	index_mapper = []
	f = open(sample_names_file)
	cell_donors = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		donor_name = line.split(':')[0]
		index_mapper.append(donor_mapping[donor_name])
		cell_donors.append(donor_name)
	f.close()
	index_mapper = np.asarray(index_mapper)
	cell_donors = np.asarray(cell_donors)
	if np.array_equal(donor_names[index_mapper],cell_donors) == False:
		print('assumption erororo')
	return index_mapper

sig_eqtl_file = sys.argv[1]
all_eqtl_file = sys.argv[2]
processed_genotype_dir = sys.argv[3]
output_root = sys.argv[4]
surge_latent_factor_file = sys.argv[5]
surge_latent_factor_num = int(sys.argv[6])
sample_names_file = sys.argv[7]


donor_genotype_to_cell_genotype_mapping = get_donor_genotype_to_cell_genotype_index_mapping(processed_genotype_dir + 'clues_immvar_chrom_' + str(1) + '.DS2.FORMAT', sample_names_file)

loadings = np.loadtxt(surge_latent_factor_file)
loading_vec = loadings[:,(surge_latent_factor_num-1)]

# Gene_name, sentinal_snp, lb, ub, savefile
test_info_file = output_root + 'test_info.txt'
t = open(test_info_file,'w')
t.write('gene_name\tchrom_num\teqtl_data_file\tld_mat_file\n')


for chrom_num in range(1,23):
	print(chrom_num)
	gene_mapping = extract_significant_genes(sig_eqtl_file, chrom_num)
	gene_mapping = map_genes_to_all_variants(all_eqtl_file, chrom_num, gene_mapping)
	variant_to_genotype_vec_mapping = make_variant_to_sample_genotype_vec_mapping(processed_genotype_dir + 'clues_immvar_chrom_' + str(chrom_num) + '.DS2.FORMAT')
	for gene_name in gene_mapping.keys():
		gene_pass_boolean, eqtl_data, ld_mat = reorganize_data_for_susie_input(gene_mapping[gene_name], variant_to_genotype_vec_mapping, loading_vec, donor_genotype_to_cell_genotype_mapping)
		if gene_pass_boolean == False:
			continue
		eqtl_data_file = output_root + gene_name + '_eqtl_data.txt'
		ld_data_file = output_root + gene_name + '_ld_mat.txt'
		t.write(gene_name + '\t' + str(chrom_num) + '\t' + eqtl_data_file + '\t' + ld_data_file + '\n')
		np.savetxt(eqtl_data_file, eqtl_data, fmt="%s", delimiter='\t')
		np.savetxt(ld_data_file, ld_mat, fmt="%s", delimiter='\t')
t.close()