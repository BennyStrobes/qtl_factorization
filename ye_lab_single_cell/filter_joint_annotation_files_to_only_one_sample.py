import numpy as np 
import os
import sys
import pdb
import gzip





def extract_sample_index(sample_names_file_parr, sample_name):
	f = open(sample_names_file_parr)
	sample_index = -1
	for line_index, line in enumerate(f):
		line_sample_name = line.rstrip()
		if line_sample_name == sample_name:
			sample_index = line_index
	f.close()
	if sample_index == -1:
		print('assumption eroror')
		pdb.set_trace()
	return sample_index



def filter_annot_file(joint_annot_file, sample_specific_annot_file, sample_index):
	f = open(joint_annot_file)
	t = open(sample_specific_annot_file,'w')
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		sample_specific_annot = 4 + sample_index
		t.write('\t'.join(data[:4]) + '\t' + data[sample_specific_annot] + '\n')
	f.close()
	t.close()

def filter_anno_ld_score_file(joint_anno_ld_score_file, sample_specific_anno_ld_score_file, sample_index):
	f = gzip.open(joint_anno_ld_score_file)
	t = gzip.open(sample_specific_anno_ld_score_file,'w')
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		sample_specific_annot = 3 + sample_index
		t.write('\t'.join(data[:3]) + '\t' + data[sample_specific_annot] + '\n')
	f.close()
	t.close()

def filter_M_file(joint_M_file, sample_specific_M_file, sample_index):
	temp = np.loadtxt(joint_M_file,dtype=str)
	t = open(sample_specific_M_file,'w')
	t.write(temp[sample_index])
	t.close()


sample_name = sys.argv[1]  # sample name
sample_names_file_parr = sys.argv[2]  # ordered list of sample names
per_sample_joint_annotation_file_stem = sys.argv[3]  # input file
per_sample_annotation_file_stem = sys.argv[4]  #output file



# Extract index corresponding to sample_namd
sample_index = extract_sample_index(sample_names_file_parr, sample_name)


# Loop through chromosomes (all annotations are per chromosome files)
for chrom_num in range(1,23):
	# Part 1: filter .annot file to just annotations from only this sample
	joint_annot_file = per_sample_joint_annotation_file_stem + '.' + str(chrom_num) + '.annot'
	sample_specific_annot_file = per_sample_annotation_file_stem + '.' + str(chrom_num) + '.annot'
	filter_annot_file(joint_annot_file, sample_specific_annot_file, sample_index)

	# Part 2: Filter .l2.ldscore.gz to just annotations from only this sample
	joint_anno_ld_score_file = per_sample_joint_annotation_file_stem + '.' + str(chrom_num) + '.l2.ldscore.gz'
	sample_specific_anno_ld_score_file = per_sample_annotation_file_stem + '.' + str(chrom_num) + '.l2.ldscore.gz'
	filter_anno_ld_score_file(joint_anno_ld_score_file, sample_specific_anno_ld_score_file, sample_index)

	# Part 3: Filter .l2.M file
	joint_M_file = per_sample_joint_annotation_file_stem + '.' + str(chrom_num) + '.l2.M'
	sample_specific_M_file = per_sample_annotation_file_stem + '.' + str(chrom_num) + '.l2.M'
	#filter_M_file(joint_M_file, sample_specific_M_file, sample_index)

	# Part 4: Filter .l2.M_5_50 file
	joint_M_file = per_sample_joint_annotation_file_stem + '.' + str(chrom_num) + '.l2.M_5_50'
	sample_specific_M_file = per_sample_annotation_file_stem + '.' + str(chrom_num) + '.l2.M_5_50'
	filter_M_file(joint_M_file, sample_specific_M_file, sample_index)

