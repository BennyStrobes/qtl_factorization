import numpy as np 
import os
import sys
import pdb
import gzip


def extract_sig_variants_from_egene_file(sig_egene_file):
	dicti = {}
	f = open(sig_egene_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_name = data[0]
		variant_name_short = variant_name.split(':')[0] + ':' + variant_name.split(':')[1]
		dicti[variant_name_short] = 1
	f.close()
	return dicti



sldsc_input_data_dir = sys.argv[1] # Directory with baseline annotations (for reference)
surge_interaction_eqtl_dir = sys.argv[2]  # directory with surge eqtls
sldsc_processed_data_dir = sys.argv[3]  # output dir


fdr = '05'
surge_eqtls = {}
for latent_factor in range(1,11):
	sig_egene_file = surge_interaction_eqtl_dir + 'surge_interaction_eqtl_cis_window_200000_factor_interaction_eqtl_results_latent_factor_' + str(latent_factor) + '_genome_wide_signficant_bf_fdr_0.' + fdr + '.txt'
	sig_variant_dicti = extract_sig_variants_from_egene_file(sig_egene_file)
	surge_eqtls[latent_factor] = sig_variant_dicti

standard_eqtl_sig_egene_file = surge_interaction_eqtl_dir + 'surge_interaction_eqtl_cis_window_200000_factor_standard_eqtl_results_genome_wide_signficant_bf_fdr_0.' + fdr + '.txt'
standard_eqtls = extract_sig_variants_from_egene_file(standard_eqtl_sig_egene_file)




for chrom_num in range(1,23):
	baseline_annot_file = sldsc_input_data_dir + 'baseline_v1.2/' + 'baseline.' + str(chrom_num) + '.annot.gz'
	surge_anno_file = sldsc_processed_data_dir + 'surge_egenes_' + fdr + '.' + str(chrom_num) + '.annot'
	joint_surge_anno_file = sldsc_processed_data_dir + 'joint_surge_egenes_' + fdr + '.' + str(chrom_num) + '.annot'
	f = gzip.open(baseline_annot_file)
	t = open(surge_anno_file, 'w')
	t_joint = open(joint_surge_anno_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			# print header
			t.write('\t'.join(data[:4]) + '\t' + 'surge_eqtl_1\tsurge_eqtl_2\tsurge_eqtl_3\tsurge_eqtl_4\tsurge_eqtl_5\tsurge_eqtl_6\tsurge_eqtl_7\tsurge_eqtl_8\tsurge_eqtl_9\tsurge_eqtl_10\tstandard_eqtl\n')
			t_joint.write('\t'.join(data[:4]) + '\t' + 'joint_surge_eqtl\tstandard_eqtl\n')
			continue
		# error checking
		if len(data) != 57:
			print('assumption eroror')
			pdb.set_trace()
		# extract relevent fields
		variant_name_short = data[0] + ':' + data[1]
		# Print header
		t.write('\t'.join(data[:4]))
		t_joint.write('\t'.join(data[:4]))
		joint_binary = '0.0'
		for latent_factor in range(1,11):
			if variant_name_short in surge_eqtls[latent_factor]:
				joint_binary = '1.0'
				t.write('\t' + '1.0')
			else:
				t.write('\t' + '0.0')
		t_joint.write('\t' + joint_binary)
		if variant_name_short in standard_eqtls:
			t.write('\t' + '1.0' + '\n')
			t_joint.write('\t' + '1.0' + '\n')
		else:
			t.write('\t' + '0.0' + '\n')
			t_joint.write('\t' + '0.0' + '\n')

	t.close()
	f.close()
	t_joint.close()




