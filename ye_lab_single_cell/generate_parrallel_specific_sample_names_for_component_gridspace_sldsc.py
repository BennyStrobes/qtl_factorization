import numpy as np 
import os
import sys
import pdb




def get_total_number_of_samples(sample_names_file):
	num_lines1 = 0
	f = open(sample_names_file)
	for line in f:
		line = line.rstrip()
		num_lines1 = num_lines1 + 1
	f.close()
	return num_lines1


# For parallelization purposes
def parallelization_start_and_end(num_tasks, job_number, total_jobs):
    tasks_per_job = (num_tasks/total_jobs) + 1
    start_task = job_number*tasks_per_job
    end_task = (job_number + 1)*tasks_per_job -1 
    return start_task, end_task


def filter_lines_in_file(input_file, output_file, start_number, end_number):
	f = open(input_file)
	t = open(output_file,'w')

	line_counter = -1
	for line_counter, line in enumerate(f):
		line = line.rstrip()
		# Skip genes not in this parallelization run
		if line_counter < start_number or line_counter > end_number:
			continue
		t.write(line + '\n')
	f.close()
	t.close()



sample_names_file = sys.argv[1]
sample_names_file_parr = sys.argv[2]
job_number = int(sys.argv[3])
num_jobs = int(sys.argv[4])


# First get total number of samples
num_samples = get_total_number_of_samples(sample_names_file)

#For parallelization purposes
start_number, end_number = parallelization_start_and_end(num_samples, job_number, num_jobs)

# Filter lines in file to [start_number, end_number] only
filter_lines_in_file(sample_names_file, sample_names_file_parr, start_number, end_number)


