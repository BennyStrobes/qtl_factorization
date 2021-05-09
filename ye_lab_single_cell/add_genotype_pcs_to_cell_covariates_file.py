import numpy as np 
import os
import sys
import pdb






cell_covariates_file = sys.argv[1]
sc_rna_seq_individual_file = sys.argv[2]
genotype_pcs_file = sys.argv[3]
processed_genotype_dir = sys.argv[4]


ordered_indis = np.loadtxt(sc_rna_seq_individual_file, dtype=str)

indi_to_genotype_pc = {}
counter = 0
f = open(genotype_pcs_file)
for line in f:
	line = line.rstrip()
	data = line.split()
	if data[0].startswith('HC'):
		line_id = data[0]
	else:
		line_id = data[0] + '_' + data[1]
	if line_id not in ordered_indis:
		print('assumption erroror')
		pdb.set_trace()
	pcs = data[2:]
	if line_id in indi_to_genotype_pc:
		print('assumptoin eororr')
		pdb.set_trace()
	indi_to_genotype_pc[line_id] = pcs
f.close()

indi_to_ancestry = {}
f = open(cell_covariates_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	pop_cov = data[4]
	ind_cov = data[1]
	indi_to_ancestry[ind_cov] = pop_cov
f.close()



####################
# Print genotype pcs to output file (individual level)
####################
individual_level_genotype_pcs = processed_genotype_dir + 'individual_level_genotype_pcs.txt'
t = open(individual_level_genotype_pcs,'w')
t.write('individual\tancestry\tgeno_pc1\tgeno_pc2\tgeno_pc3\tgeno_pc4\n')
for indi in ordered_indis:
	t.write(indi + '\t' + indi_to_ancestry[indi] + '\t' + '\t'.join(indi_to_genotype_pc[indi][:4]) + '\n')
t.close()

####################
# Print genotype pcs to output file (cell level)
####################
cell_level_genotype_pcs = processed_genotype_dir + 'pseudobulk_sample_covariates_with_genotype_pcs.txt'
t = open(cell_level_genotype_pcs,'w')


f = open(cell_covariates_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\tgenotype_pc1\tgenotype_pc2\n')
		continue
	indi_id = data[1]
	genotype_pcs = indi_to_genotype_pc[indi_id]
	t.write(line + '\t' + genotype_pcs[0] + '\t' + genotype_pcs[1] + '\n')
f.close()
t.close()
