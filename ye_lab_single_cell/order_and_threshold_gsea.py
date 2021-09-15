import numpy as np 
import os
import sys
import pdb









####################
# Command line args
####################
input_file = sys.argv[1]


# Parse output file from input file
output_file = input_file.split('.')[0] + '_significant_sorted.txt'



f = open(input_file)
t = open(output_file,'w')

count = 0
head_count = 0
pathways = []
for line in f:
	line = line.rstrip()
	if count < 3:
		# skip first 3 lines
		count = count + 1
		continue
	if head_count == 0:
		# Print header
		head_count = head_count + 1
		data = line.split('\t')
		t.write('\t'.join(data[1:]) + '\n')
		continue
	data = line.split('\t')
	bf_corrected_pval = float(data[7])
	pathways.append((bf_corrected_pval,line))
f.close()
# Sort
sorted_pathways = sorted(pathways, key=lambda tup: tup[0])
for pathway in sorted_pathways:
	if pathway[0] < .05:
		t.write(pathway[1] + '\n')
t.close()