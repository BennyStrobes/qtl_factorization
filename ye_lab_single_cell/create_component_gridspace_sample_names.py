import numpy as np 
import os
import sys
import pdb






loading_file = sys.argv[1]
sample_names_file = sys.argv[2]
num_samples = int(sys.argv[3])

t = open(sample_names_file,'w')

loadings = np.loadtxt(loading_file)

num_components = loadings.shape[1]


for component_num in range(num_components):
	miny = np.min(loadings[:, component_num])
	maxy = np.max(loadings[:, component_num])
	size = maxy - miny

	gridspace = np.arange(miny, maxy, size/num_samples)

	for position in gridspace:
		t.write('component_' + str(component_num) + ':' + str(position) + '\n')

t.close()