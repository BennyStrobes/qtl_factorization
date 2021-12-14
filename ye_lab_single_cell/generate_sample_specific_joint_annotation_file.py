import numpy as np 
import os
import sys
import pdb
import gzip






def create_mapping_from_variant_to_per_sample_effects(variant_level_sample_specific_eqtl_effect_sizes_file):
	f = open(variant_level_sample_specific_eqtl_effect_sizes_file)
	dicti = {}
	head_count = 0
	counter = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			sample_names = np.asarray(data[1:])
			continue
		variant_name = data[0]
		betas = np.asarray(data[1:])
		variant_name_short = variant_name.split(':')[0] + ':' + variant_name.split(':')[1]
		# quick erroro checking
		if variant_name_short in dicti:
			counter = counter + 1
			continue
		dicti[variant_name_short] = betas
	f.close()

	return dicti, sample_names








variant_level_sample_specific_eqtl_effect_sizes_file = sys.argv[1]
sldsc_input_data_dir = sys.argv[2]
per_sample_joint_annotation_file_stem = sys.argv[3]

# First create mapping from variant names to array of per sample squared effects sizes
variant_to_per_sample_effects, sample_names = create_mapping_from_variant_to_per_sample_effects(variant_level_sample_specific_eqtl_effect_sizes_file)


zero_annotation_string = '\t'.join(np.zeros(len(sample_names)).astype(str))



for chrom_num in range(1,23):
	print(chrom_num)
	baseline_annot_file = sldsc_input_data_dir + 'baseline_v1.2/' + 'baseline.' + str(chrom_num) + '.annot.gz'
	per_sample_anno_file = per_sample_joint_annotation_file_stem + '.' + str(chrom_num) + '.annot'
	f = gzip.open(baseline_annot_file)
	t = open(per_sample_anno_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			# print header
			t.write('\t'.join(data[:4]) + '\t' + '\t'.join(sample_names) + '\n')
			continue
		# error checking
		if len(data) != 57:
			print('assumption eroror')
			pdb.set_trace()
		# extract relevent fields
		variant_name_short = data[0] + ':' + data[1]
		# Print header
		t.write('\t'.join(data[:4]) + '\t')
		if variant_name_short not in variant_to_per_sample_effects:
			t.write(zero_annotation_string + '\n')
		else:
			t.write('\t'.join(variant_to_per_sample_effects[variant_name_short]) + '\n')
	t.close()
	f.close()
