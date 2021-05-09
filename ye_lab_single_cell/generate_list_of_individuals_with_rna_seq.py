import numpy as np 
import os
import sys
import pdb




###################
# Command line args
###################
cell_covariates_file = sys.argv[1]
sc_rna_seq_individual_file = sys.argv[2]

# Extract dictionary set of indiviudals from covariate file
indis = {}
f = open(cell_covariates_file)
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		continue
	indis[data[1]] = 1
f.close()

# Print individuals to output file
t = open(sc_rna_seq_individual_file, 'w')
for indi in sorted(indis.keys()):
	t.write(indi + '\n')
t.close()