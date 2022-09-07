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

def map_genes_to_all_variants(all_eqtl_file, chrom_num, gene_mapping, eqtl_data_version):
	f = open(all_eqtl_file)
	chrom_string = str(chrom_num) + ':'
	variant_to_gwas = {}
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
		if eqtl_data_version == 'version1':
			beta = data[2]
			std_err_beta = data[3]
			var_beta = str(np.square(float(std_err_beta)))
			pvalue = data[4]
		elif eqtl_data_version == 'version2':
			beta = data[3]
			std_err_beta = data[4]
			var_beta = str(np.square(float(std_err_beta)))
			pvalue = data[2]
		else:
			print('assumption erororor')
		gene_mapping[gene_name].append((variant_id, beta, var_beta, pvalue))
		variant_to_gwas[variant_id] = (-1, -1, -1)
		variant_to_gwas[variant_id_alt] = (-1, -1, -1)
	f.close()
	return gene_mapping, variant_to_gwas

def map_variants_to_gwas_summary_stats(gwas_file, variant_to_gwas):
	f = open(gwas_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 5:
			print('assumption erroor')
			pdb.set_trace()
		variant_id = data[0]
		if variant_id not in variant_to_gwas:
			continue
		variant_to_gwas[variant_id] = (float(data[1]), float(data[2]), float(data[3]))
	f.close()
	return variant_to_gwas

def merge_study_pair(eqtl_variants, variant_to_gwas):
	variant_names = []
	snp_positions = []
	eqtl_betas = []
	eqtl_beta_vars = []
	eqtl_pvalues = []
	gwas_betas = []
	gwas_beta_vars = []
	gwas_pvalues = []
	for eqtl_variant_tuple in eqtl_variants:
		eqtl_variant = eqtl_variant_tuple[0]
		variant_info = eqtl_variant.split(':')
		eqtl_variant_alt = variant_info[0] + ':' + variant_info[1] + ':' + variant_info[3] + ':' + variant_info[2]
		eqtl_variant_beta = float(eqtl_variant_tuple[1])
		eqtl_variant_beta_var = float(eqtl_variant_tuple[2])
		eqtl_variant_pvalue = float(eqtl_variant_tuple[3])
		if variant_to_gwas[eqtl_variant][1] != -1 and variant_to_gwas[eqtl_variant_alt][1] != -1:
			print('assumption error')
			pdb.set_trace()
		elif variant_to_gwas[eqtl_variant][1] != -1 or variant_to_gwas[eqtl_variant_alt][1] != -1:
			if variant_to_gwas[eqtl_variant][1] != -1:
				gwas_beta = variant_to_gwas[eqtl_variant][0]
				gwas_beta_var = variant_to_gwas[eqtl_variant][1]
				gwas_pvalue = variant_to_gwas[eqtl_variant][2]
			elif variant_to_gwas[eqtl_variant_alt][1] != -1:
				gwas_beta = variant_to_gwas[eqtl_variant_alt][0]
				gwas_beta_var = variant_to_gwas[eqtl_variant_alt][1]
				gwas_pvalue = variant_to_gwas[eqtl_variant_alt][2]
			else:
				print('assumption eroror')
				pdb.set_trace()
			variant_names.append(eqtl_variant)
			snp_positions.append(eqtl_variant.split(':')[1])
			eqtl_betas.append(eqtl_variant_beta)
			eqtl_beta_vars.append(eqtl_variant_beta_var)
			eqtl_pvalues.append(eqtl_variant_pvalue)
			gwas_betas.append(gwas_beta)
			gwas_beta_vars.append(gwas_beta_var)
			gwas_pvalues.append(gwas_pvalue)
	if len(gwas_pvalues) == 0:
		return False, np.nan, np.nan
	elif np.min(np.asarray(gwas_pvalues)) < 1e-3 and len(gwas_pvalues) > 50:
		# Put into organized matrix (coloc format)
		eqtl_mat = np.transpose(np.vstack((np.asarray(variant_names), np.asarray(snp_positions), np.asarray(eqtl_betas).astype(str), np.asarray(eqtl_beta_vars).astype(str), np.asarray(eqtl_pvalues).astype(str))))
		gwas_mat = np.transpose(np.vstack((np.asarray(variant_names), np.asarray(snp_positions), np.asarray(gwas_betas).astype(str), np.asarray(gwas_beta_vars).astype(str), np.asarray(gwas_pvalues).astype(str))))
		# Add headers
		eqtl_mat_w_header = np.vstack((['snp', 'position', 'beta', 'varbeta', 'pvalue'], eqtl_mat))
		gwas_mat_w_header = np.vstack((['snp', 'position', 'beta', 'varbeta', 'pvalue'], gwas_mat))
		return True, eqtl_mat_w_header, gwas_mat_w_header
	else:
		return False, np.nan, np.nan


sig_eqtl_file = sys.argv[1]
all_eqtl_file = sys.argv[2]
gwas_study_file_root = sys.argv[3]
output_root = sys.argv[4]
eqtl_data_version = sys.argv[5]



# Gene_name, sentinal_snp, lb, ub, savefile
test_info_file = output_root + 'test_info.txt'
t = open(test_info_file,'w')
t.write('gene_name\tchrom_num\teqtl_data_file\tgwas_data_file\n')


for chrom_num in range(1,23):
	print(chrom_num)
	gene_mapping = extract_significant_genes(sig_eqtl_file, chrom_num)
	gene_mapping, variant_to_gwas = map_genes_to_all_variants(all_eqtl_file, chrom_num, gene_mapping, eqtl_data_version)
	pdb.set_trace()
	variant_to_gwas = map_variants_to_gwas_summary_stats(gwas_study_file_root + str(chrom_num) + '.txt', variant_to_gwas)

	for gene_name in gene_mapping.keys():
		gene_pass_boolean, eqtl_data, gwas_data = merge_study_pair(gene_mapping[gene_name], variant_to_gwas)
		if gene_pass_boolean == False:
			continue
		eqtl_data_file = output_root + gene_name + '_eqtl_data.txt'
		gwas_data_file = output_root + gene_name + '_gwas_data.txt'
		t.write(gene_name + '\t' + str(chrom_num) + '\t' + eqtl_data_file + '\t' + gwas_data_file + '\n')
		np.savetxt(eqtl_data_file, eqtl_data, fmt="%s", delimiter='\t')
		np.savetxt(gwas_data_file, gwas_data, fmt="%s", delimiter='\t')
t.close()
