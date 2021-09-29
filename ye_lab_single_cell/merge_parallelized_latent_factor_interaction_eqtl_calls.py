import numpy as np 
import os
import sys
import pdb
import scipy.stats

def bf_fdr_multiple_testing_correction(variant_gene_pairs_eqtl_results_file, multple_testing_correction_results_file, fdr_thresh):
	f = open(variant_gene_pairs_eqtl_results_file)
	t = open(multple_testing_correction_results_file, 'w')
	t.write('snp_id\tgene_id\tbeta\tstd_err_beta\tpvalue\tnum_snps_in_gene\tfdr\n')
	head_count = 0
	genes = {}
	for line in f:
		line = line.rstrip()
		data = line.split()
		gene_id = data[1]
		variant_id = data[0]
		pvalue = float(data[4])
		beta = float(data[2])
		std_err = float(data[3])
		abs_t_value = np.abs(beta/std_err)
		if gene_id not in genes:
			genes[gene_id] = (variant_id, pvalue, abs_t_value, 1, line)
		else:
			old_t_value = genes[gene_id][2]
			old_count = genes[gene_id][3]
			if abs_t_value >= old_t_value:
				genes[gene_id] = (variant_id, pvalue, abs_t_value, old_count+1, line)
			else:
				genes[gene_id] = (genes[gene_id][0], genes[gene_id][1], genes[gene_id][2], old_count+1, genes[gene_id][4])
	f.close()
	# Loop through genes and do BF correction
	bf_gene_array = []
	for gene in genes.keys():
		lead_variant = genes[gene][0]
		lead_nominal_pvalue = genes[gene][1]
		lead_nominal_abs_tvalue = genes[gene][2]
		num_variants_at_gene = genes[gene][3]
		test_line = genes[gene][4]
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
			line = gene_tuple[4]
			data = line.split('\t')
			std_err_beta = (data[3])
			t.write(data[0] + '\t' + data[1] + '\t' + data[2]  + '\t' + std_err_beta + '\t' + data[4] + '\t' + str(gene_tuple[3]) + '\t' + str(fdr) + '\n')
	t.close()

def make_sure_files_exist(output_root, total_jobs, suffix):
	booly = True
	for job_number in range(total_jobs):
		file_name = output_root + str(job_number) + '_' + str(total_jobs) + '_results' + suffix
		if os.path.isfile(file_name) == False:
			print(file_name)
			booly = False
	return booly

def merge_parallelized_results(output_root, suffix, total_jobs, lf_num):
	to_run = make_sure_files_exist(output_root, total_jobs, suffix)
	if to_run == False:
		print('Missing required input files. Please re-evaluate :)')
		return
	# Open output (merged result) file handle
	t = open(output_root + 'latent_factor_' + str(lf_num+1) + '_merged' + suffix, 'w')
	t2 = open(output_root + 'latent_factor_' + str(lf_num+1) + '_merged_include_nan' + suffix, 'w')
	# Loop through parrallelized jobs
	for job_number in range(total_jobs):
		file_name = output_root + str(job_number) + '_' + str(total_jobs) + '_results' + suffix
		# Open file for one job
		f = open(file_name)
		# To identify header
		head_count = 0
		# Stream file from one job
		counter = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			beta = float(data[(2 + lf_num*3)])
			std_err_beta = float(data[(3 + lf_num*3)])
			old_pvalue = data[(4 + lf_num*3)]
			new_pvalue = scipy.stats.norm.cdf(-np.abs( beta/std_err_beta))*2
			new_line = data[0] + '\t' + data[1] + '\t' + str(beta) + '\t' + str(std_err_beta) + '\t' + str(new_pvalue)
			new_data = new_line.split('\t')
			t2.write(new_line + '\n')
			counter = counter +1
			if len(data) < (lf_num*3):
				print('miss')
				continue
			if new_data[2] == 'NA':
				continue
			# HEADER
			#if head_count == 0:
			#	head_count = head_count + 1
				# Print header if this the first job
			#	if job_number == 0:
			#		t.write(line + '\n')
			#	continue
			# Standard line
			t.write(new_line + '\n')
		f.close()
		# Delete file from single job
		#os.system ('rm ' + file_name)
	t2.close()
	t.close()

def get_number_of_latent_factors(fdr_file):
	f = open(fdr_file)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			continue
		num_lf = len(data[3].split(','))
		break
	f.close()
	return num_lf

def number_of_hits_per_latent_factor(fdr_file, num_hits_per_lf_file, nominal_coefficient_pvalue_thresholds):
	# open output file handle
	t = open(num_hits_per_lf_file, 'w')
	# Print header
	t.write('latent_factor\tnominal_pvalue_threshold\tnumber_of_hits\n')
	# Extract number of latent factors
	num_latent_factors = get_number_of_latent_factors(fdr_file)
	# Looop through nominal p-value thresholds
	for nominal_pvalue in nominal_coefficient_pvalue_thresholds:
		# Initialize counter to keep track of number of hits per latent factors
		num_hits = np.zeros(num_latent_factors)
		f = open(fdr_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			coef_pvalz = np.asarray(data[3].split(',')).astype(float)
			num_hits = num_hits + 1.0*(coef_pvalz < nominal_pvalue)
		f.close()
		for latent_factor_num in range(num_latent_factors):
			t.write(str(latent_factor_num+1) + '\t' + str(nominal_pvalue) + '\t' + str(num_hits[latent_factor_num]) + '\n')
	t.close()

def get_quantile_of_pvalue(pvalue, quantile_values):
	num_quantiles = len(quantile_values) -1

	pvalue_quantile = -2
	hits = 0
	for quantile_num in range(num_quantiles):
		if pvalue >= quantile_values[quantile_num] and pvalue < quantile_values[(quantile_num+1)]:
			pvalue_quantile = quantile_num
			hits = hits +1
	if hits == 0 and pvalue == quantile_values[-1]:
		pvalue_quantile = len(quantile_values) - 3
		hits = hits + 1
	if hits != 1:
		print('assumption error')
		pdb.set_trace()
	return pvalue_quantile

def number_of_hits_per_latent_factor_seperated_by_quantiles(fdr_file, num_hits_per_lf_file, num_quantiles):
	temp = np.loadtxt(fdr_file, dtype=str, delimiter='\t')
	pvalues = temp[1:,2].astype(float)

	quantile_size = 1.0/num_quantiles
	quantiles = np.arange(0.0, 1.00000000001, quantile_size)
	quantile_values = []
	for quantile in quantiles:
		quantile_values.append(np.quantile(pvalues, quantile))
	quantile_values = np.asarray(quantile_values)

	nominal_pvalue = 1e-4

	# open output file handle
	t = open(num_hits_per_lf_file, 'w')
	# Print header
	t.write('latent_factor\tquantile_bin\tnumber_of_hits\n')
	# Extract number of latent factors
	num_latent_factors = get_number_of_latent_factors(fdr_file)
	# Looop through nominal p-value thresholds
	for quantile_num in range(num_quantiles):
		# Initialize counter to keep track of number of hits per latent factors
		num_hits = np.zeros(num_latent_factors)
		f = open(fdr_file)
		head_count = 0
		for line in f:
			line = line.rstrip()
			data = line.split('\t')
			if head_count == 0:
				head_count = head_count + 1
				continue
			pvalue = float(data[2])
			pvalue_quantile = get_quantile_of_pvalue(pvalue, quantile_values)
			if pvalue_quantile != quantile_num:
				continue
			coef_pvalz = np.asarray(data[3].split(',')).astype(float)
			num_hits = num_hits + 1.0*(coef_pvalz < nominal_pvalue)
		f.close()
		for latent_factor_num in range(num_latent_factors):
			t.write(str(latent_factor_num+1) + '\t' + str(quantile_num) + '\t' + str(num_hits[latent_factor_num]) + '\n')
	t.close()

def get_number_of_latent_factors(file_name):
	f = open(file_name)
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split('\t')
		if head_count == 0:
			head_count = head_count + 1
			num_lf = int((len(data)-2.0)/3)
			continue
	f.close()
	return num_lf

output_root = sys.argv[1]
total_jobs = int(sys.argv[2])

# Extract number of latent factors
num_lf = get_number_of_latent_factors(output_root + '0_' + str(total_jobs) + '_results.txt')


for lf_num in range(num_lf):
	print(lf_num)
	merged_file = output_root + 'latent_factor_' + str(lf_num+1) + '_merged.txt'
	merge_parallelized_results(output_root, '.txt', total_jobs, lf_num)

	fdr = .05
	fdr_file = output_root + 'latent_factor_' + str(lf_num+1) + '_genome_wide_signficant_bf_fdr_' + str(fdr) + '.txt'
	bf_fdr_multiple_testing_correction(merged_file, fdr_file, fdr)

	fdr = .1
	fdr_file = output_root + 'latent_factor_' + str(lf_num+1) + '_genome_wide_signficant_bf_fdr_' + str(fdr) + '.txt'
	bf_fdr_multiple_testing_correction(merged_file, fdr_file, fdr)

	fdr = .2
	fdr_file = output_root + 'latent_factor_' + str(lf_num+1) + '_genome_wide_signficant_bf_fdr_' + str(fdr) + '.txt'
	bf_fdr_multiple_testing_correction(merged_file, fdr_file, fdr)

	fdr = .5
	fdr_file = output_root + 'latent_factor_' + str(lf_num+1) + '_genome_wide_signficant_bf_fdr_' + str(fdr) + '.txt'
	bf_fdr_multiple_testing_correction(merged_file, fdr_file, fdr)


