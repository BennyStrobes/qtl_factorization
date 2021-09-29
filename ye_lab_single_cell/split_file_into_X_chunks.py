import numpy as np 
import os
import sys
import pdb





input_file = sys.argv[1]
total_lines = int(sys.argv[2]) - 1
num_chunks = int(sys.argv[4])
output_stem = sys.argv[5]
header = sys.argv[6]


lines_per_job = int(np.ceil(total_lines/num_chunks))

job_counter = 0
line_counter = 0
head_count = 0

t = open(output_stem + str(job_counter) + '.txt', 'w')

f = open(input_file)
for line in f:
	if head_count == 0 and header == "True":
		head_count = head_count + 1
		continue
	t.write(line)
	line_counter = line_counter + 1
	if np.mod(line_counter, lines_per_job) == 0:
		t.close()
		job_counter = job_counter + 1
		t = open(output_stem + str(job_counter) + '.txt', 'w')
t.close()
f.close()