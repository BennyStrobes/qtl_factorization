import numpy as np 
import os
import sys
import pdb






variant_gene_level_eqtl_effect_sizes = sys.argv[1]
variant_level_eqtl_effect_sizes = sys.argv[2]


# Mapping from variant to vector of efffect sizes across samples
dicti = {}
f = open(variant_gene_level_eqtl_effect_sizes)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		sample_ids = np.asarray(data[2:])
		num_samples = len(sample_ids)
		continue
	variant_id = data[0]
	betas = np.asarray(data[2:]).astype(float)
	if variant_id not in dicti:
		dicti[variant_id] = np.square(betas)
	else:
		dicti[variant_id] = dicti[variant_id] + np.square(betas)
f.close()


# Stream to output file
t = open(variant_level_eqtl_effect_sizes,'w')
t.write('variant_id\t' + '\t'.join(sample_ids) + '\n')

for variant_id in sorted(dicti.keys()):
	squared_summed_betas = dicti[variant_id]
	t.write(variant_id + '\t' + '\t'.join(squared_summed_betas.astype(str)) + '\n')
t.close()