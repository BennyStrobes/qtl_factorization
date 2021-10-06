import numpy as np 
import os
import sys
import pdb




def extract_egene_variant_gene_pairs(egene_file, num_lf, variant_gene_pairs, standard_eqtl_variant_gene_pairs):
	f = open(egene_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_name = data[0]
		gene_name = data[1]
		test_name = variant_name + '_' + gene_name
		variant_gene_pairs[test_name] = np.zeros(num_lf) - 1.0
		standard_eqtl_variant_gene_pairs[test_name] = -1.0
	f.close()
	return variant_gene_pairs, standard_eqtl_variant_gene_pairs

def fill_in_pvalues_to_variant_gene_pair_dictionary(eqtl_file, lf_num, variant_gene_pairs):
	f = open(eqtl_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		variant_name = data[0]
		gene_name = data[1]
		test_name = variant_name + '_' + gene_name
		if test_name not in variant_gene_pairs:
			continue
		#pvalue = float(data[4])
		beta = float(data[2])
		std_err = float(data[3])
		abs_t = np.abs(beta/std_err)
		variant_gene_pairs[test_name][(lf_num-1)] = abs_t
	f.close()
	return variant_gene_pairs

def fill_in_pvalues_to_standard_eqtl_variant_gene_pair_dictionary(eqtl_file, s_variant_gene_pairs):
	f = open(eqtl_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		variant_name = data[0]
		gene_name = data[1]
		test_name = variant_name + '_' + gene_name
		if test_name not in s_variant_gene_pairs:
			continue
		#pvalue = float(data[4])
		beta = float(data[2])
		std_err = float(data[3])
		abs_t = np.abs(beta/std_err)
		s_variant_gene_pairs[test_name] = abs_t
	f.close()
	return s_variant_gene_pairs

surge_interaction_eqtl_dir = sys.argv[1]
eqtl_correlation_dir = sys.argv[2]


U_file = surge_interaction_eqtl_dir + 'surge_interaction_eqtl_cis_window_200000_factor_surge_latent_factors.txt'
U = np.loadtxt(U_file)
loading_corr = np.abs(np.corrcoef(np.transpose(U)))
np.savetxt(eqtl_correlation_dir + 'loading_correlation.txt', loading_corr, fmt="%s", delimiter='\t')


# First extract list of variant gene pairs
# Our list will be list of variant-gene pairs that are egenes
variant_gene_pairs = {}
standard_eqtl_variant_gene_pairs = {}
for lf_num in range(1,11):
	egene_file = surge_interaction_eqtl_dir + 'surge_interaction_eqtl_cis_window_200000_factor_interaction_eqtl_results_latent_factor_' + str(lf_num) + '_genome_wide_signficant_bf_fdr_0.05.txt'
	variant_gene_pairs, standard_eqtl_variant_gene_pairs = extract_egene_variant_gene_pairs(egene_file, 10, variant_gene_pairs, standard_eqtl_variant_gene_pairs)

# Now fill in dictionary of variant-gene pairs with p-values
for lf_num in range(1,11):
	eqtl_file = surge_interaction_eqtl_dir + 'surge_interaction_eqtl_cis_window_200000_factor_interaction_eqtl_results_latent_factor_' + str(lf_num) + '_merged.txt'
	variant_gene_pairs = fill_in_pvalues_to_variant_gene_pair_dictionary(eqtl_file, lf_num, variant_gene_pairs)

eqtl_file = surge_interaction_eqtl_dir + 'surge_interaction_eqtl_cis_window_200000_factor_standard_eqtl_results_merged.txt'
standard_eqtl_variant_gene_pairs = fill_in_pvalues_to_standard_eqtl_variant_gene_pair_dictionary(eqtl_file, standard_eqtl_variant_gene_pairs)

test_names = np.asarray([*variant_gene_pairs])
num_tests = len(test_names)

sum_stat_mat = np.zeros((num_tests, 10))

for i, test_name in enumerate(test_names):
	test_arr = variant_gene_pairs[test_name]
	# simple error checking
	if sum(test_arr == -1) != 0:
		print('assumption errror')
		pdb.set_trace()
	sum_stat_mat[i,:] = test_arr

cor_mat = np.corrcoef(np.transpose(sum_stat_mat))
np.savetxt(eqtl_correlation_dir + 'eqtl_correlation.txt', cor_mat, fmt="%s", delimiter='\t')


sum_stat_mat = np.zeros((num_tests, 11))

for i, test_name in enumerate(test_names):
	test_arr = variant_gene_pairs[test_name]
	test_arr = np.hstack((test_arr,[standard_eqtl_variant_gene_pairs[test_name]]))
	# simple error checking
	if sum(test_arr == -1) != 0:
		print('assumption errror')
		pdb.set_trace()
	sum_stat_mat[i,:] = test_arr


cor_mat = np.corrcoef(np.transpose(sum_stat_mat))
np.savetxt(eqtl_correlation_dir + 'all_eqtl_correlation.txt', cor_mat, fmt="%s", delimiter='\t')






