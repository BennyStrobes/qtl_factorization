import numpy as np 
import os
import sys
import pdb

def bf_fdr_multiple_testing_correction(variant_gene_pairs_eqtl_results_file, multple_testing_correction_results_file, fdr_thresh):
	f = open(variant_gene_pairs_eqtl_results_file)
	t = open(multple_testing_correction_results_file, 'w')
	head_count = 0
	genes = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\tnum_snps_in_gene\tfdr\n')
			continue
		gene_id = data[1]
		variant_id = data[0]
		if data[7] == 'NA':
			data[7] = '1.0'
		pvalue = float(data[7])

		if gene_id not in genes:
			genes[gene_id] = (variant_id, pvalue, 1, line)
		else:
			old_pvalue = genes[gene_id][1]
			old_count = genes[gene_id][2]
			if pvalue <= old_pvalue:
				genes[gene_id] = (variant_id, pvalue, old_count+1, line)
			else:
				genes[gene_id] = (genes[gene_id][0], genes[gene_id][1], old_count+1, genes[gene_id][3])
	f.close()
	# Loop through genes and do BF correction
	bf_gene_array = []
	for gene in genes.keys():
		lead_variant = genes[gene][0]
		lead_nominal_pvalue = genes[gene][1]
		num_variants_at_gene = genes[gene][2]
		test_line = genes[gene][3]
		bf_corrected_pvalue = lead_nominal_pvalue*num_variants_at_gene
		if bf_corrected_pvalue > 1.0:
			bf_corrected_pvalue = 1.0
		bf_gene_array.append((bf_corrected_pvalue, lead_variant, gene, num_variants_at_gene, test_line))
	sorted_bf_gene_array = sorted(bf_gene_array, key=lambda tup: tup[0])
	# BH correction
	kk = 1
	num_genes = len(sorted_bf_gene_array)
	sig = True
	for gene_tuple in sorted_bf_gene_array:
		bf_pvalue = gene_tuple[0]
		fdr = num_genes*bf_pvalue/kk 
		kk = kk + 1
		if fdr > fdr_thresh:
			sig = False
		if sig == True:
			t.write(gene_tuple[4] + '\t' + str(gene_tuple[3]) + '\t' + str(fdr) + '\n')
	t.close()

def make_sure_files_exist(output_root, total_jobs, suffix):
	booly = True
	for job_number in range(total_jobs):
		file_name = output_root + str(job_number) + '_' + str(total_jobs) + '_results' + suffix
		if os.path.isfile(file_name) == False:
			print(file_name)
			booly = False
	return booly

def merge_parallelized_results(output_root, suffix, total_jobs):
	to_run = make_sure_files_exist(output_root, total_jobs, suffix)
	if to_run == False:
		print('Missing required input files. Please re-evaluate :)')
		return
	# Open output (merged result) file handle
	t = open(output_root + 'merged' + suffix, 'w')
	t.write('variant_id\tgene_id\tpvalue\tbeta\n')
	# Loop through parrallelized jobs
	for job_number in range(total_jobs):
		file_name = output_root + str(job_number) + '_' + str(total_jobs) + '_results' + suffix
		# Open file for one job
		f = open(file_name)
		# To identify header
		head_count = 0
		# Stream file from one job
		for line in f:
			line = line.rstrip()
			# HEADER
			#if head_count == 0:
			#	head_count = head_count + 1
				# Print header if this the first job
			#	if job_number == 0:
			#		t.write(line + '\n')
			#	continue
			# Standard line
			data = line.split('\t')
			#if data[2] == 'NA':
				#continue
			t.write(line + '\n')
		f.close()
		# Delete file from single job
		#os.system('rm ' + file_name)
	t.close()

def get_pvalue_only_file(input_file, output_file):
	f = open(input_file)
	t = open(output_file,'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		if data[7] == 'NA':
			data[7] = '1.0'
		t.write(data[7] + '\n')
	f.close()
	t.close()


output_root = sys.argv[1]
total_jobs = int(sys.argv[2])


merged_file = output_root + 'merged.txt'
merge_parallelized_results(output_root, ".txt", total_jobs)



merged_file = output_root + 'merged.txt'
merged_pvalue_only_file = output_root + 'merged_pvalue_only.txt'
get_pvalue_only_file(merged_file,merged_pvalue_only_file)

fdr = .05
fdr_file = output_root + 'genome_wide_signficant_bf_fdr_' + str(fdr) + '.txt'
bf_fdr_multiple_testing_correction(merged_file, fdr_file, fdr)
