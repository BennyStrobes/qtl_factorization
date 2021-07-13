import numpy as np 
import os
import sys
import pdb


def get_frq_file_variant_names(frq_file):
	f = open(frq_file)
	used_variants = {}
	head_count = 0
	arr1 = []
	arr2 = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if len(data) != 6:
			print('assumtpoineori eororo')
			pdb.set_trace()
		variant_name1 = data[0] + ':' + data[1]
		a1_info = data[4].split(':')
		a2_info = data[5].split(':')
		if len(a1_info) != 2:
			print('assumption erorro')
			pdb.set_trace()
		if len(a2_info) != 2:
			print('assumption error')
			pdb.set_trace()
		ref_allele = a1_info[0]
		alt_allele = a2_info[0]
		variant_name2 = data[0] + ':' + data[1] + ':' + ref_allele + ':' + alt_allele
		if variant_name2 in used_variants:
			print('assumption error')
			pdb.set_trace()
		used_variants[variant_name2] = 1
		arr1.append(variant_name1)
		arr2.append(variant_name2)
	return np.asarray(arr1), np.asarray(arr2)


#######################
# Command line args
#######################
frq_file = sys.argv[1]
dosage_file = sys.argv[2]
dosage2_file = sys.argv[3]


variant_names1, variant_names2 = get_frq_file_variant_names(frq_file)

f = open(dosage_file)
t = open(dosage2_file, 'w')

if len(variant_names2) != len(np.unique(variant_names2)):
	print('assumption erororor')
	pdb.set_trace()

head_count = 0
counter = 0
for line in f:
	line = line.rstrip()
	data = line.split('\t')
	if head_count == 0:
		head_count = head_count + 1
		num_header_features = len(data)
		t.write(data[0] + '\t' + data[1] + '\t' + 'VARIANT_ID\t' + '\t'.join(data[2:]) + '\n')
		continue
	if len(data) != num_header_features:
		print('assumption eororro')
	line_variant_id = data[0] + ':' + data[1]
	if line_variant_id != variant_names1[counter]:
		print('assumption erororo')
		pdb.set_trace()
	t.write(data[0] + '\t' + data[1] + '\t' + variant_names2[counter] + '\t' + '\t'.join(data[2:]) + '\n')
	counter = counter + 1
f.close()
t.close()
