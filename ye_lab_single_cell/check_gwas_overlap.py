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



def get_gwas_study_info(gwas_files):
	f = open(gwas_files)
	head_count = 0
	counter = 0
	study_names = []
	study_paths = []
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		study_name = data[0]
		study_path = data[1]
		study_names.append(study_name)
		study_paths.append(study_path)
	f.close()
	return np.asarray(study_names), np.asarray(study_paths)

def extract_eqtls(egene_file):
	f = open(egene_file)
	head_count = 0
	egenes = {}
	evariants = {}
	eqtls = {}
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		snp_id = data[0]
		gene_id = data[1]
		pvalue = data[4]
		egenes[gene_id] = 1
		evariants[snp_id] = 1
		eqtls[gene_id + '_' + snp_id] = pvalue
	f.close()
	return egenes, evariants, eqtls

def create_rs_snp_mappings(gwas_snp_info_dir):
	snp_to_rs = {}
	for chrom_num in range(1,23):
		snp_file = gwas_snp_info_dir + '1000G.EUR.QC.' + str(chrom_num) + '.bim'
		f = open(snp_file)
		for line in f:
			line = line.rstrip()
			data = line.split()
			if len(data) != 6:
				print('assumption error')
				pdb.set_trace()
			snp_id_1 = data[0] + ':' + data[3] + ':' + data[4] + ':' + data[5]
			snp_id_2 = data[0] + ':' + data[3] + ':' + data[5] + ':' + data[4]
			rs_id = data[1]
			if snp_id_1 in snp_to_rs or snp_id_2 in snp_to_rs:
				print('assumption eroror')
				pdb.set_trace()
			snp_to_rs[snp_id_1] = rs_id
			snp_to_rs[snp_id_2] = rs_id
		f.close()
	return snp_to_rs

def get_pvalues_of_overlapping_variants(rs_id_to_pvalue_mapping, rs_id_to_snp_id):
	pvalues = []
	for rs_id in rs_id_to_snp_id.keys():
		if rs_id in rs_id_to_pvalue_mapping:
			pvalues.append(rs_id_to_pvalue_mapping[rs_id])
	return np.asarray(pvalues)

# Return the bin number corresponding to this distance
def get_distance_bin(distance, distance_bin_size):
    return int(math.floor(distance/distance_bin_size))


# Return the bin number corresponding to this distance
def get_maf_bin(maf, maf_bin_size):
    return int(math.floor(maf/maf_bin_size))

def make_background_object(test_info_file):
	distance_bin_size = 10000
	maf_bin_size = .05
	eqtl_distance = 50000
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
	return background_variants


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


def get_rsids_from_snp_ids(evariants, snp_to_rs):
	rs_ids = {}
	for evariant in evariants.keys():
		if evariant in snp_to_rs:
			rs_ids[snp_to_rs[evariant]] = 1
	return rs_ids

def get_study_rs_id_to_pvalue_mapping(gwas_study_file):
	dicti = {}
	f = open(gwas_study_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		rs_id = data[0]
		if len(data) == 6:
			pvalue =chi2.sf(np.square(float(data[5])), 1)
		if len(data) == 4:
			pvalue =chi2.sf(np.square(float(data[3])), 1)
		if rs_id in dicti:
			print('assumptione oeoror')
			pdb.set_trace()
		dicti[rs_id] = pvalue
	f.close()
	return dicti

def get_study_snp_id_to_pvalue_mapping(gwas_study_file_root):
	dicti = {}
	for chrom_num in range(1,23):
		file_name = gwas_study_file_root + str(chrom_num) + '.txt'
		f = open(file_name)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			snp_id = data[0]
			snp_id_info = snp_id.split(':')
			snp_id_alt = snp_id_info[0] + ':' + snp_id_info[1] + ':' + snp_id_info[3] + ':' + snp_id_info[2]
			pvalue = float(data[3])
			if snp_id in dicti or snp_id_alt in dicti:
				print('assumption eororor')
			dicti[snp_id] = pvalue
			dicti[snp_id_alt] = pvalue
		f.close()
	return dicti

egene_file = sys.argv[1]
gwas_files = sys.argv[2]
test_info_file = sys.argv[3]
gwas_snp_info_dir = sys.argv[4]
output_root = sys.argv[5]




# Extract names and files of gwas studies
gwas_study_names, gwas_study_files = get_gwas_study_info(gwas_files)

immune_indices = np.asarray([3, 5, 7, 8, 9, 10, 11, 14])
non_immune_indices = np.asarray([0,1,2, 4, 6, 12, 13])


# Create output files
output_file = output_root + 'cross_study_gwas_enrichment.txt'
t = open(output_file,'w')
t.write('\t'.join(gwas_study_names) + '\n')

gwas_snp_id_to_pvalue = []
for i, gwas_study_name in enumerate(gwas_study_names):
	print(i)
	gwas_study_file_root = gwas_study_files[i]
	gwas_snp_id_to_pvalue.append(get_study_snp_id_to_pvalue_mapping(gwas_study_file_root))
#f = open('mypickle.pickle', 'wb')
#pickle.dump(gwas_snp_id_to_pvalue, f)
#f.close()

#f = open('mypickle.pickle', 'rb')
#gwas_snp_id_to_pvalue = pickle.load(f)

# Extract egenes and evariants for tests
egenes, evariants, eqtls = extract_eqtls(egene_file)


test_pi1s = []
for i, gwas_study_name in enumerate(gwas_study_names):
	pvalues = get_pvalues_of_overlapping_variants(gwas_snp_id_to_pvalue[i], evariants)
	pi0 = qvalue.estimate(pvalues)
	pi1 = 1.0 - pi0
	test_pi1s.append(pi1)
test_pi1s = np.asarray(test_pi1s)
#test_pi1s_immune_average = np.mean(test_pi1s[immune_indices])
#test_pi1s_non_immune_average = np.mean(test_pi1s[non_immune_indices])




######################
# Simulate background variants matched for maf and distTSS
num_background_runs = 1000

bgrd_object = make_background_object(test_info_file)

bgrd_variants = get_bgrd_eqtls(num_background_runs, bgrd_object, eqtls, test_info_file)

for background_run in range(num_background_runs):
	print(background_run)

	bgrd_run_evariants = bgrd_variants[background_run]
	# Get test rs-ids
	#bgrd_rs_ids = get_rsids_from_snp_ids(bgrd_run_evariants, snp_to_rs)

	bgrd_pi1s = []

	for i, gwas_study_name in enumerate(gwas_study_names):
		bgrd_pvalues = get_pvalues_of_overlapping_variants(gwas_snp_id_to_pvalue[i], bgrd_run_evariants)
		pi0 = qvalue.estimate(bgrd_pvalues)
		pi1 = 1.0 - pi0
		bgrd_pi1s.append(pi1)

	bgrd_pi1s = np.asarray(bgrd_pi1s)
	#bgrd_pi1s_immune_average = np.mean(bgrd_pi1s[immune_indices])
	#bgrd_pi1s_non_immune_average = np.mean(bgrd_pi1s[non_immune_indices])
	delta_pi1s = test_pi1s - bgrd_pi1s
	#print(ranksum_stats)
	t.write('\t'.join((delta_pi1s).astype(str)) + '\n')
	t.flush()
t.close()
