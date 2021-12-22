import numpy as np 
import os
import sys
import pdb








eqtl_effect_sizes_input_file = sys.argv[1]
eqtl_effect_sizes_output_file = sys.argv[2]




# Open input file handle
f = open(eqtl_effect_sizes_input_file)
# Open output file handle
t = open(eqtl_effect_sizes_output_file,'w')
# Write header to output file
t.write('variant_id\tgene_id\tstatic_eqtl\n')

# Stream input file (each line is a variant-gene pair)
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	# Extract relevent fields
	variant_id = data[0]
	gene_id = data[1]

	beta_shared = data[2]

	# Print to output file
	t.write(variant_id + '\t' + gene_id + '\t' + beta_shared + '\n')

f.close()
t.close()