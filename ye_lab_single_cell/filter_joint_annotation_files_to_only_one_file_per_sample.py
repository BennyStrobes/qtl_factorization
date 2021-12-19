import numpy as np 
import os
import sys
import pdb
import gzip






def get_sample_names(sample_names_file):
	f = open(sample_names_file)
	arr = []
	for line in f:
		line = line.rstrip()
		arr.append(line)
	f.close()
	return np.asarray(arr)


def generate_annot_files_for_each_sample(per_sample_joint_annotation_file_stem, sample_names):
	arr = []  # each line in arr is a snp
	for chrom_num in range(1,23):
		chrom_annotation_file = per_sample_joint_annotation_file_stem + '.' + str(chrom_num) + '.annot'
		f = open(chrom_annotation_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				# Quick error checking
				if np.array_equal(np.asarray(data[4:]), sample_names) == False:
					print('assumption eroror')
					pdb.set_trace()
				continue
			annot_across_samples = np.asarray(data[4:]).astype(float)
			arr.append(annot_across_samples)
		f.close()

	# Get into matrix where each row is a snp and each column is a sample
	arr = np.array(arr)

	# QUick error check
	if arr.shape[1] != len(sample_names):
		print('assumption erroror')
		pdb.set_trace()

	for sample_iter, sample_name in enumerate(sample_names):
		# Get annotation specific to this sample
		sample_specific_annotation = arr[:, sample_iter]
		# Name of file to save this annotation
		sample_specific_annotation_npy_file = per_sample_joint_annotation_file_stem + '_samplename_' + sample_name + '_annot.npy'
		np.save(sample_specific_annotation_npy_file, sample_specific_annotation)

def generate_annot_ld_score_files_for_each_sample(per_sample_joint_annotation_file_stem, sample_names):
	arr = []  # each line in arr is a snp
	for chrom_num in range(1,23):
		chrom_annotation_file = per_sample_joint_annotation_file_stem + '.' + str(chrom_num) + '.l2.ldscore.gz'
		f = gzip.open(chrom_annotation_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				# Quick error checking
				un_processed_header = np.asarray(data[3:])
				processed_header = []
				for ele in un_processed_header:
					processed_header.append(ele.split('L2')[0])
				if np.array_equal(np.asarray(processed_header), sample_names) == False:
					print('assumption eroror')
					pdb.set_trace()
				continue
			annot_across_samples = np.asarray(data[3:]).astype(float)
			arr.append(annot_across_samples)
		f.close()

	# Get into matrix where each row is a snp and each column is a sample
	arr = np.array(arr)

	# QUick error check
	if arr.shape[1] != len(sample_names):
		print('assumption erroror')
		pdb.set_trace()

	for sample_iter, sample_name in enumerate(sample_names):
		# Get annotation specific to this sample
		sample_specific_annotation = arr[:, sample_iter]
		# Name of file to save this annotation
		sample_specific_annotation_npy_file = per_sample_joint_annotation_file_stem + '_samplename_' + sample_name + '_ldscore.npy'
		np.save(sample_specific_annotation_npy_file, sample_specific_annotation)

def generate_M_files_for_each_sample(per_sample_joint_annotation_file_stem, sample_names, suffix):
	arr = []  # each line in arr is a snp
	for chrom_num in range(1,23):
		chrom_M_file = per_sample_joint_annotation_file_stem + '.' + str(chrom_num) + '.l2.' + suffix
		data = np.loadtxt(chrom_M_file)
		arr.append(data)
	arr = np.asarray(arr)
	# Take sum across chromosomes
	per_sample_M = np.sum(arr,axis=0)

	for sample_iter, sample_name in enumerate(sample_names):
		# Get annotation specific to this sample
		sample_M = per_sample_M[sample_iter]
		# Name of file to save this annotation
		sample_specific_annotation_npy_file = per_sample_joint_annotation_file_stem + '_samplename_' + sample_name + '_l2_' + suffix + '.npy'
		np.save(sample_specific_annotation_npy_file, sample_M)


sample_names_file = sys.argv[1]
per_sample_joint_annotation_file_stem = sys.argv[2]
per_sample_annotation_summary_file = sys.argv[3]


sample_names = get_sample_names(sample_names_file)


# Generate annot files for each sample
generate_annot_files_for_each_sample(per_sample_joint_annotation_file_stem, sample_names)

# Generate annot ld scorefiles for each sample
generate_annot_ld_score_files_for_each_sample(per_sample_joint_annotation_file_stem, sample_names)

generate_M_files_for_each_sample(per_sample_joint_annotation_file_stem, sample_names, 'M')
generate_M_files_for_each_sample(per_sample_joint_annotation_file_stem, sample_names, 'M_5_50')

# Print to summary file
t = open(per_sample_annotation_summary_file,'w')
for sample_name in sample_names:
	t.write(sample_name)
	t.write('\t' + per_sample_joint_annotation_file_stem + '_samplename_' + sample_name + '_annot.npy')
	t.write('\t' + per_sample_joint_annotation_file_stem + '_samplename_' + sample_name + '_ldscore.npy')
	t.write('\t' + per_sample_joint_annotation_file_stem + '_samplename_' + sample_name + '_l2_' + 'M' + '.npy')
	t.write('\t' + per_sample_joint_annotation_file_stem + '_samplename_' + sample_name + '_l2_' + 'M_5_50' + '.npy')
	t.write('\n')
t.close()



