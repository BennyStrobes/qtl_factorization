import numpy as np 
import os
import sys
import pdb







#####################
# Command line args
#####################
xcell_sample_enrichment_file = sys.argv[1]
cell_type = sys.argv[2]
cell_type_context_file = sys.argv[3]

# Load in sample data
sample_data = np.loadtxt(xcell_sample_enrichment_file,dtype=str,delimiter='\t')
sample_data_header = sample_data[0,:]
column_indices = np.where(sample_data_header==cell_type)[0]

if len(column_indices) != 1:
	print('assumption eorororor')
	pdb.set_trace()

column_index = column_indices[0]

xcell_enrichments = sample_data[1:,column_index]

np.savetxt(cell_type_context_file, xcell_enrichments, fmt="%s", delimiter='\t')