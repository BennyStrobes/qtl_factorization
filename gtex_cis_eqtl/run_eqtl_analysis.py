import numpy as np 
import os
import sys
import pdb
import statsmodels.api as sm
from sklearn import linear_model



# Generate eqtl effect size (beta) and p-value for one test
def run_eqtl_one_test(expression, genotype, covariates):
	# Covariate matrix
	X = np.vstack((genotype,covariates.T)).T
	# Add intercept
	X2 = sm.add_constant(X)
	est = sm.OLS(expression, X2)
	est2 = est.fit()
	beta = est2.params[1]
	standard_error = est2.bse[1]
	pvalue = est2.pvalues[1]

	return beta, standard_error, pvalue


# Run eQTL analysis
def eqtl_analysis(covariate_file, test_names_file, expression_file, genotype_file, variant_gene_pairs_eqtl_results_file):
	# Load in covariates (fixed across all tests)
	covariates = np.loadtxt(covariate_file)
	# Open up file handles
	test_name_handle = open(test_names_file)
	expression_handle = open(expression_file)
	genotype_handle = open(genotype_file)
	# Output file handle
	t = open(variant_gene_pairs_eqtl_results_file, 'w')

	# Loop through tests
	head_count = 0  # Used to skip header
	counter = 0
	for line in test_name_handle:
		test_name_line  = line.rstrip()
		# Skip header
		if head_count == 0:
			head_count = head_count + 1
			t.write(test_name_line + '\tbeta\tstandard_error\tpvalue\n')
			continue
		counter = counter + 1
		if np.mod(counter, 1000) == 0.0:
			print(counter)
		expression = np.asarray(expression_handle.next().rstrip().split('\t')).astype(float)
		genotype = np.asarray(genotype_handle.next().rstrip().split('\t')).astype(float)
		beta, std_err, pvalue = run_eqtl_one_test(expression, genotype, covariates)
		t.write(test_name_line + '\t' + str(beta) + '\t' + str(std_err) + '\t' + str(pvalue) + '\n')
	t.close()
	test_name_handle.close()
	expression_handle.close()
	genotype_handle.close()

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
		gene_id = data[0]
		variant_id = data[1]
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

def bf_top_nn_multiple_testing_correction(variant_gene_pairs_eqtl_results_file, multple_testing_correction_results_file, nn_thresh):
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
		gene_id = data[0]
		variant_id = data[1]
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
		if kk > (nn_thresh+1):
			sig = False
		if sig == True:
			t.write(gene_tuple[4] + '\t' + str(gene_tuple[3]) + '\t' + str(fdr) + '\n')
	t.close()



#####################
# Command line args
#####################
expression_file = sys.argv[1]  # Input file of dimension num_testsXnum_samples
genotype_file = sys.argv[2]  # Input file of dimension num_testsXnum_samples
test_names_file = sys.argv[3]
covariate_file = sys.argv[4]  # Input file of dimension num_samplesXnum_covariates
output_root = sys.argv[5]  # output root


####################
# Run eQTL analysis
####################
# Output file
variant_gene_pairs_eqtl_results_file = output_root + 'all_variant_gene_pairs.txt'
eqtl_analysis(covariate_file, test_names_file, expression_file, genotype_file, variant_gene_pairs_eqtl_results_file)


####################
# Multiple-testing correction
####################
# Output file
fdr_thresh=.01
multple_testing_correction_results_file = output_root + 'multiple_testing_bf_bh_' + str(fdr_thresh) + '_fdr_' + '.txt'
# Perform bonferonni correction at variant level (per gene) and then BH at gene level
#bf_fdr_multiple_testing_correction(variant_gene_pairs_eqtl_results_file, multple_testing_correction_results_file, fdr_thresh)


####################
# Multiple-testing correction
####################
# Output file
fdr_thresh=.05
multple_testing_correction_results_file = output_root + 'multiple_testing_bf_bh_' + str(fdr_thresh) + '_fdr_' + '.txt'
# Perform bonferonni correction at variant level (per gene) and then BH at gene level
#bf_fdr_multiple_testing_correction(variant_gene_pairs_eqtl_results_file, multple_testing_correction_results_file, fdr_thresh)


####################
# Multiple-testing correction
####################
# Output file
fdr_thresh=.1
multple_testing_correction_results_file = output_root + 'multiple_testing_bf_bh_' + str(fdr_thresh) + '_fdr_' + '.txt'
# Perform bonferonni correction at variant level (per gene) and then BH at gene level
#bf_fdr_multiple_testing_correction(variant_gene_pairs_eqtl_results_file, multple_testing_correction_results_file, fdr_thresh)


####################
# Top nn genes
####################
# Output file
nn_thresh= 800
multple_testing_correction_results_file = output_root + 'multiple_testing_bf_top_' + str(nn_thresh) + '_genes.txt'
# Perform bonferonni correction at variant level (per gene) and then BH at gene level
#bf_top_nn_multiple_testing_correction(variant_gene_pairs_eqtl_results_file, multple_testing_correction_results_file, nn_thresh)