import numpy as np 
import os 
import sys
import pdb
import gzip 




input_file = sys.argv[1]
output_file = sys.argv[2]


f = gzip.open(input_file)
t = gzip.open(output_file,'w')
head_count = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		t.write(line + '\n')
		continue
	pvalue = float(data[9])
	if pvalue < 1e-300:
		continue
	t.write(line + '\n')
f.close()
t.close()