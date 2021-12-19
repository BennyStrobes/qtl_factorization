import numpy as np 
import os
import sys
import pdb
import gzip



def filter_annotation_file(input_annot_file, output_annot_file, num_annot, annot_to_keep):
	f = gzip.open(input_annot_file)
	t = gzip.open(output_annot_file,'w')
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		header = np.asarray(data[:4])
		annot = np.asarray(data[4:])
		if len(annot) != num_annot:
			print('assumption eroor')
			pdb.set_trace()
		t.write('\t'.join(header) + '\t' + '\t'.join(annot[annot_to_keep]) + '\n')
	f.close()
	t.close()

def filter_ldscore_file(input_annot_file, output_annot_file, num_annot, annot_to_keep):
	f = gzip.open(input_annot_file)
	t = gzip.open(output_annot_file,'w')
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		header = np.asarray(data[:3])
		annot = np.asarray(data[3:])
		if len(annot) != num_annot:
			print('assumption eroor')
			pdb.set_trace()
		t.write('\t'.join(header) + '\t' + '\t'.join(annot[annot_to_keep]) + '\n')
	f.close()
	t.close()


def filter_m_file(input_m_file, output_m_file, num_annot, annot_to_keep):
	data = np.loadtxt(input_m_file)
	if len(data) != num_annot:
		print('assumption eroror')
		pdb.set_trace()
	t = open(output_m_file,'w')
	new_data = data[annot_to_keep]
	t.write('\t'.join(new_data.astype(str)) + '\n')
	t.close()



input_stem = sys.argv[1]
output_stem = sys.argv[2]

num_annot = 97
annot_to_remove = np.asarray([72, 73, 74, 75])
annot_to_keep = []
for annot_num in range(num_annot):
	if annot_num not in annot_to_remove:
		annot_to_keep.append(annot_num)
annot_to_keep = np.asarray(annot_to_keep)

for chrom_num in range(1,23):
	print(chrom_num)
	# Annotation file
	input_annot_file = input_stem + str(chrom_num) + '.annot.gz'
	output_annot_file = output_stem + str(chrom_num) + '.annot.gz'
	filter_annotation_file(input_annot_file, output_annot_file, num_annot, annot_to_keep)
	
	# LD score file
	input_ldscore_file = input_stem + str(chrom_num) + '.l2.ldscore.gz'
	output_ldscore_file = output_stem + str(chrom_num) + '.l2.ldscore.gz'
	filter_ldscore_file(input_ldscore_file, output_ldscore_file, num_annot, annot_to_keep)

	# M file
	input_m_file = input_stem + str(chrom_num) + '.l2.M'
	output_m_file = output_stem + str(chrom_num) + '.l2.M'
	filter_m_file(input_m_file, output_m_file, num_annot, annot_to_keep)

	# M_5_50 file
	input_m_file = input_stem + str(chrom_num) + '.l2.M_5_50'
	output_m_file = output_stem + str(chrom_num) + '.l2.M_5_50'
	filter_m_file(input_m_file, output_m_file, num_annot, annot_to_keep)

