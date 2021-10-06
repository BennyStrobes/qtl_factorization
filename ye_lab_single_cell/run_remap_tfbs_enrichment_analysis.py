import numpy as np 
import os
import sys
import pdb
from scipy.stats import chi2
import math
import random
import scipy.stats
import pickle
import qvalue

def make_chromosome_with_remap_tfbs(chrom_num, remap_tfbs_file):
	chrom_num_string = 'chr' + str(chrom_num) + '\t'
	chromosome = ['NULL']*259250621
	f = open(remap_tfbs_file)
	counter = 0
	for line in f:
		# Ignore peaks not on current chromosome
		if line.startswith(chrom_num_string) == False:
			continue
		line = line.rstrip()
		data = line.split()
		# Error checking
		if len(data) != 9:
			print('assumption erororor')
			pdb.set_trace()
		# Extract relevent fields
		start = int(data[1])
		end = int(data[2])
		peak_name = data[3].split(':')[0]
		# More error checking
		if end < start:
			print('assumption error')
			pdb.set_trace()
		for pos in range(start, end):
			if chromosome[pos] == 'NULL':
				chromosome[pos] = peak_name
			else:
				chromosome[pos] = chromosome[pos] + ';' + peak_name
	f.close()
	return chromosome

def parse_peak_name(peak_name):
	arr = []
	peak_info = peak_name.split(':')
	# error checking
	if len(peak_info) != 2:
		print('assumption error')
		pdb.set_trace()
	gene_name = peak_info[0]
	cell_types = peak_info[1].split(',')
	for cell_type in cell_types:
		arr.append(gene_name + ':' + cell_type)
	return arr

def get_tf_dicti(remap_tfbs_file):
	tf_dicti = {}
	chrom_num_string = 'chr' + str(1) + '\t'
	f = open(remap_tfbs_file)
	for line in f:
		# Ignore peaks not on current chromosome
		#if line.startswith(chrom_num_string) == False:
			#continue
		line = line.rstrip()
		data = line.split()
		# Extract relevent fields
		start = int(data[1])
		end = int(data[2])
		peak_name = data[3]
		info = peak_name.split(':')
		if len(info) != 2:
			print('assumption eroror')
			pdb.set_trace()
		tf_name = info[0]
		tf_dicti[tf_name] = 0
		#unique_tfs = parse_peak_name(peak_name)
		#for tf_name in unique_tfs:
			#tf_dicti[tf_name] = 0
	f.close()
	return tf_dicti

def get_sig_variants(surge_interaction_sig_eqtl_file):
	f = open(surge_interaction_sig_eqtl_file)
	variants = []
	eqtls = {}
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		variant_id = data[0]
		gene_id = data[1]
		variants.append(variant_id)
		eqtls[gene_id + '_' + variant_id] = 1
	f.close()
	return np.unique(np.asarray(variants)), eqtls

def update_tf_dicti(tf_dicti, sig_variants, chromosome, chrom_num):
	chrom_string = str(chrom_num)
	for variant in sig_variants:
		line_chrom_num = variant.split(':')[0]
		if line_chrom_num != chrom_string:
			continue
		pos = int(variant.split(':')[1])
		info = chromosome[pos]
		if info == 'NULL':
			continue
		infos = info.split(';')
		for peak_name in infos:
			tf_dicti[peak_name] = tf_dicti[peak_name] + 1
			'''
			unique_tfs = parse_peak_name(peak_name)
			for unique_tf in unique_tfs:
				tf_dicti[unique_tf] = tf_dicti[unique_tf] + 1
			'''
	return tf_dicti

def update_bgrd_tf_dicti(bgrd_tf_dicti, bgrd_variants, bgrd_run_num, chromosome, chrom_num):
	chrom_string = str(chrom_num)
	for variant in bgrd_variants:
		line_chrom_num = variant.split(':')[0]
		if line_chrom_num != chrom_string:
			continue
		pos = int(variant.split(':')[1])
		info = chromosome[pos]
		if info == 'NULL':
			continue
		infos = info.split(';')
		for peak_name in infos:
			bgrd_tf_dicti[peak_name][bgrd_run_num] = bgrd_tf_dicti[peak_name][bgrd_run_num] + 1
			'''
			unique_tfs = parse_peak_name(peak_name)
			for unique_tf in unique_tfs:
				bgrd_tf_dicti[unique_tf][bgrd_run_num] = bgrd_tf_dicti[unique_tf][bgrd_run_num] + 1
			'''
	return bgrd_tf_dicti

# Return the bin number corresponding to this distance
def get_distance_bin(distance, distance_bin_size):
    return int(math.floor(distance/distance_bin_size))


# Return the bin number corresponding to this distance
def get_maf_bin(maf, maf_bin_size):
    return int(math.floor(maf/maf_bin_size))

def make_background_object(test_info_file):
	distance_bin_size = 10000
	maf_bin_size = .05
	eqtl_distance = 200000
	####################
	# Initialize object
	####################
	background_qtls = []
	# number of bins needed for maf and distance
	num_distance_bins = int(math.ceil(eqtl_distance/distance_bin_size + 1))
	num_maf_bins = int(math.ceil(.5/maf_bin_size + 1))
	# Add each possible bin
	for distance_bin in range(num_distance_bins):
		background_qtls.append([])
		for maf_bin in range(num_maf_bins):
			background_qtls[distance_bin].append([])
	f = open(test_info_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		snp_id = data[1]
		test_name = gene_id + '_' + snp_id
		distance = float(data[5])
		maf = float(data[6])
		# Return the bin number corresponding to this distance
		distance_bin = get_distance_bin(distance, distance_bin_size)
		# Return the bin number corresponding to this distance
		maf_bin = get_maf_bin(maf, maf_bin_size)
		background_qtls[distance_bin][maf_bin].append(test_name)
	f.close()
	return background_qtls

def sample_background_variant_gene_pairs(test_infos, bgrd_object):
	distance_bin_size = 10000
	maf_bin_size = .05
	background_variant_gene_pairs = {}
	for test_info in test_infos:
		test_distance = test_info[0]
		test_maf = test_info[1]
		distance_bin = get_distance_bin(test_distance, distance_bin_size)
		maf_bin = get_maf_bin(test_maf, maf_bin_size)
		converged = False
		while converged == False:
			if len(bgrd_object[distance_bin][maf_bin]) < 4:
				print('small backgroudn: should investigate')
				print(len(bgrd_object[distance_bin][maf_bin]))
			randomly_selected_pair = random.choice(bgrd_object[distance_bin][maf_bin])
			if randomly_selected_pair not in background_variant_gene_pairs:
				background_variant_gene_pairs[randomly_selected_pair] = 1
				converged = True
	if len(background_variant_gene_pairs) != len(test_infos):
		print('assumption eorororro')
		pdb.set_trace()
	background_variants = {}
	for test in background_variant_gene_pairs.keys():
		variant_id = test.split('_')[1]
		background_variants[variant_id] = 1
	return np.asarray([*background_variants])


def get_bgrd_eqtls(num_background_runs, bgrd_object, eqtls, test_info_file):
	distance_bin_size = 10000
	maf_bin_size = .05

	# Extract ordered list of maf and distance of each tests
	test_infos = []
	f = open(test_info_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			continue
		gene_id = data[0]
		snp_id = data[1]
		test_name = gene_id + '_' + snp_id
		distance = float(data[5])
		maf = float(data[6])
		if test_name not in eqtls:
			continue
		# This test is an eqtl
		# Return the bin number corresponding to this distance
		distance_bin = get_distance_bin(distance, distance_bin_size)
		# Return the bin number corresponding to this distance
		maf_bin = get_maf_bin(maf, maf_bin_size)
		test_infos.append((distance, maf))
		#bin_tests = bgrd_object[distance_bin][maf_bin]
	f.close()
	background_variant_gene_pairs = []
	for perm_num in range(num_background_runs):
		background_variant_gene_pairs.append(sample_background_variant_gene_pairs(test_infos, bgrd_object))
	return background_variant_gene_pairs

def extract_test_and_bgrd_genes_for_gsea(output_file, output_enrichment_file, output_probability_file, gsea_test_genes_file, gsea_bgrd_genes_file, num_background_runs):
	f = open(output_file)
	t = open(output_enrichment_file,'w')
	t.write('TF\todds_ratio\n')
	arr = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		tf = data[0]
		aa = float(data[1])
		bb = float(data[2]) 
		cc = np.asarray(data[3].split(',')).astype(float) 
		dd = np.asarray(data[4].split(',')).astype(float) 
		orat = ((aa+1)/(bb+1))/((cc+1)/(dd+1))
		for ele in orat:
			t.write(tf + '\t' + str(ele) + '\n')
		count = sum(aa < cc)
		if aa == 0:
			arr.append((data[0], count, 0.0))
		else:
			if np.isnan(np.nanmean(aa/cc)):
				pdb.set_trace()
			meany = np.mean(np.ma.masked_invalid(aa/cc))
			arr.append((data[0], len(cc)-count, meany))
	f.close()
	t.close()
	t = open(output_probability_file, 'w')
	t.write('TF\tpvalue\n')
	for ele in arr:
		tf = ele[0]
		p_val = (num_background_runs - ele[1])/num_background_runs
		t.write(tf + '\t' + str(p_val) + '\n')
	t.close()
	#arr.sort(key=lambda x: x[2])
	arr.sort(key=lambda x: x[1])

	t = open(gsea_test_genes_file,'w')
	test_count = 0
	print(arr[((-201)):-1])
	for tupler in arr[((-201)):-1]:
		test_count = test_count + 1
		t.write(tupler[0] + '\n')
	t.close()
	print(str(test_count) + ' test genes')
	t = open(gsea_bgrd_genes_file,'w')
	bgrd_count = 0
	for tupler in arr[:500]:
		t.write(tupler[0] + '\n')
		bgrd_count = bgrd_count + 1
	t.close()	
	print(str(bgrd_count) + ' background genes')

surge_interaction_sig_eqtl_file = sys.argv[1]
remap_tfbs_file = sys.argv[2]
output_file = sys.argv[3]
output_enrichment_file = sys.argv[4]
output_probability_file = sys.argv[5]
test_info_file = sys.argv[6]
gsea_test_genes_file = sys.argv[7]
gsea_bgrd_genes_file = sys.argv[8]

num_background_runs = 10000
'''
sig_variants, eqtls = get_sig_variants(surge_interaction_sig_eqtl_file)

bgrd_object = make_background_object(test_info_file)

bgrd_variants = get_bgrd_eqtls(num_background_runs, bgrd_object, eqtls, test_info_file)

num_background_variants = []
for bgrd_variant_arr in bgrd_variants:
	num_background_variants.append(len(bgrd_variant_arr))
num_background_variants = np.asarray(num_background_variants)


tf_dicti = get_tf_dicti(remap_tfbs_file)
bgrd_tf_dicti = {}
for tf_name in sorted(tf_dicti.keys()):
	bgrd_tf_dicti[tf_name] = np.zeros(num_background_runs)


t = open(output_file, 'w')
for chrom_num in range(1,23):
	print(chrom_num)
	chromosome = make_chromosome_with_remap_tfbs(chrom_num, remap_tfbs_file)
	tf_dicti = update_tf_dicti(tf_dicti, sig_variants, chromosome, chrom_num)
	for bgrd_run_num in range(num_background_runs):
		bgrd_tf_dicti = update_bgrd_tf_dicti(bgrd_tf_dicti, bgrd_variants[bgrd_run_num], bgrd_run_num, chromosome, chrom_num)
for tf_name in sorted(tf_dicti.keys()):
	real_counts = tf_dicti[tf_name]
	num_real_variants = len(sig_variants)
	bgrd_counts = bgrd_tf_dicti[tf_name]
	t.write(tf_name + '\t' + str(real_counts) + '\t' + str(num_real_variants) + '\t' + ','.join(bgrd_counts.astype(str)) + '\t' + ','.join(num_background_variants.astype(str)) + '\n') 
t.close()
'''
extract_test_and_bgrd_genes_for_gsea(output_file, output_enrichment_file, output_probability_file, gsea_test_genes_file, gsea_bgrd_genes_file, num_background_runs)





