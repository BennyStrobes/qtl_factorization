import numpy as np 
import os
import sys
import pdb



def extract_enrichment_and_tau_from_bottom_line_of_sldsc_results_file(per_cell_trait_file):
	aa = np.loadtxt(per_cell_trait_file,dtype=str, delimiter='\t')
	return aa[-1,4], aa[-1,-1]




sample_names_file = sys.argv[1]
per_cell_sldsc_results_dir = sys.argv[2]

trait_names = ['ukbb_blood_monocyte_count', 'Ulcerative_Colitis']

output_file = per_cell_sldsc_results_dir + 'per_cell_sldsc_results.txt'

t = open(output_file,'w')
t.write('sample_name')
for trait_name in trait_names:
	t.write('\t' + trait_name + '_enrichment\t' + trait_name + '_tau_star')
t.write('\n')


f = open(sample_names_file)
for line in f:
	sample_name = line.rstrip()
	t.write(sample_name)
	for trait_name in trait_names:
		per_cell_trait_file = per_cell_sldsc_results_dir + 'sample_specific_eqtl_effect_sizes_' + trait_name + '_' + sample_name + '.results'
		if os.path.exists(per_cell_trait_file) == False:
			t.write('\tNA\tNA')
		else:
			enrichment, tau = extract_enrichment_and_tau_from_bottom_line_of_sldsc_results_file(per_cell_trait_file)
			t.write('\t' + enrichment + '\t' + tau)
	t.write('\n')

f.close()
t.close()
