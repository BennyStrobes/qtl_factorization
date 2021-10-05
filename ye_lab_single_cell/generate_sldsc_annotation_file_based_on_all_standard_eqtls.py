import numpy as np 
import os
import sys
import pdb
import gzip


def extract_sig_variants_from_eqtl_file(sig_egene_file, pvalue_thresh):
	dicti = {}
	f = open(sig_egene_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		pvalue = float(data[4])
		if pvalue > pvalue_thresh:
			continue
		variant_name = data[0]
		variant_name_short = variant_name.split(':')[0] + ':' + variant_name.split(':')[1]
		dicti[variant_name_short] = 1
	f.close()
	return dicti



sldsc_input_data_dir = sys.argv[1] # Directory with baseline annotations (for reference)
surge_interaction_eqtl_dir = sys.argv[2]  # directory with surge eqtls
sldsc_processed_data_dir = sys.argv[3]  # output dir


pvalue_thresh = 1e-5
surge_eqtls = {}
for latent_factor in range(1,11):
	all_eqtl_file = surge_interaction_eqtl_dir + 'surge_interaction_eqtl_cis_window_200000_factor_interaction_eqtl_results_latent_factor_' + str(latent_factor) + '_merged.txt'
	sig_variant_dicti = extract_sig_variants_from_eqtl_file(all_eqtl_file, pvalue_thresh)
	surge_eqtls[latent_factor] = sig_variant_dicti

standard_eqtl_file = surge_interaction_eqtl_dir + 'surge_interaction_eqtl_cis_window_200000_factor_standard_eqtl_results_merged.txt'
standard_eqtls = extract_sig_variants_from_eqtl_file(standard_eqtl_file, pvalue_thresh)



for chrom_num in range(1,23):
	baseline_annot_file = sldsc_input_data_dir + 'baseline_v1.2/' + 'baseline.' + str(chrom_num) + '.annot.gz'
	surge_anno_file = sldsc_processed_data_dir + 'standard_eqtls_' + str(pvalue_thresh) + '.' + str(chrom_num) + '.annot'
	f = gzip.open(baseline_annot_file)
	t = open(surge_anno_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			# print header
			t.write('\t'.join(data[:4]) + '\t' + 'standard_eqtl\n')
			continue
		# error checking
		if len(data) != 57:
			print('assumption eroror')
			pdb.set_trace()
		# extract relevent fields
		variant_name_short = data[0] + ':' + data[1]
		# Print header
		t.write('\t'.join(data[:4]))
		if variant_name_short in standard_eqtls:
			t.write('\t' + '1.0' + '\n')
		else:
			t.write('\t' + '0.0' + '\n')

	t.close()
	f.close()




