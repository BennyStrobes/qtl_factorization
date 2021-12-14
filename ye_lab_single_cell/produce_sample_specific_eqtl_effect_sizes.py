import numpy as np 
import os
import sys
import pdb








sample_names_file = sys.argv[1]
loading_file = sys.argv[2]
surge_eqtl_effect_sizes_file = sys.argv[3]
sample_specific_eqtl_effect_sizes_file = sys.argv[4]



# Load in sample names
sample_names = np.loadtxt(sample_names_file, dtype=str)

# Load in SURGE Loadings
loadings = np.loadtxt(loading_file, delimiter='\t')

# number of samples
num_samples = len(sample_names)

# Number of factors
num_factors = loadings.shape[1]

# Quick error checking
if num_samples != loadings.shape[0]:
	print('assumption eoororor')
	pdb.set_trace()


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
	# Quick error checking
	if len(betas) != num_factors:
		print('assumption erroro')
		pdb.set_trace()
	# Now compute eqtl effect sizes in each sample
	beta_in_each_sample = np.dot(loadings,betas) + beta_shared

	# Print to output file
	t.write(variant_id + '\t' + gene_id + '\t' + '\t'.join(beta_in_each_sample.astype(str)) + '\n')

f.close()
t.close()