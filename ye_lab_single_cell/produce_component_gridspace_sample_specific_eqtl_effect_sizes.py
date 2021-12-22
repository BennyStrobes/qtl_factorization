import numpy as np 
import os
import sys
import pdb








sample_names_file = sys.argv[1]
surge_eqtl_effect_sizes_file = sys.argv[2]
sample_specific_eqtl_effect_sizes_file = sys.argv[3]


pvalue_thresh=1e-6

# Load in sample names
sample_names = np.loadtxt(sample_names_file, dtype=str)

# Put sample names into information form
sample_info_arr = []
for sample_name in sample_names:
	component_num = int(sample_name.split(':')[0].split('_')[1])
	component_position = float(sample_name.split(':')[1])
	sample_info_arr.append((component_num, component_position))

# Quick error checking
if len(sample_info_arr) != len(sample_names):
	print('assumptoino eroror')
	pdb.set_trace()


# number of samples
num_samples = len(sample_names)



# Open input file handle
f = open(surge_eqtl_effect_sizes_file)
# Open output file handle
t = open(sample_specific_eqtl_effect_sizes_file,'w')
# Write header to output file
t.write('variant_id\tgene_id\t' + '\t'.join(sample_names) + '\n')

# Stream input file (each line is a variant-gene pair)
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	# Extract relevent fields
	variant_id = data[0]
	gene_id = data[1]

	beta_shared = float(data[2])
	betas = np.asarray(data[3:]).astype(float)


	beta_in_each_sample = []
	for sample_num in range(num_samples):
		sample_info = sample_info_arr[sample_num]
		sample_component_num = sample_info[0]
		sample_component_position = sample_info[1]
		sample_beta = (beta_shared + (betas[sample_component_num]*sample_component_position))
		beta_in_each_sample.append(sample_beta)
	beta_in_each_sample = np.asarray(beta_in_each_sample)

	# Print to output file
	t.write(variant_id + '\t' + gene_id + '\t' + '\t'.join(beta_in_each_sample.astype(str)) + '\n')

f.close()
t.close()