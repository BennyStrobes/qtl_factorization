import numpy as np 
import os
import sys
import pdb
import surge.surge_inference
from statsmodels.stats import diagnostic



# Load in sample overlap data
def load_in_sample_overlap_data(sample_overlap_file):
	f = open(sample_overlap_file)
	Z1 = []
	for line in f:
		line = line.rstrip()
		Z1.append(int(line))
	f.close()
	num_individuals = len(np.unique(Z1))
	return np.asarray(Z1), int(num_individuals)

def standardize_variance_ratio_between_expression_and_genotype(Y, G):
	num_tests = Y.shape[1]
	for test_num in range(num_tests):
		test_sdev = np.std(Y[:, test_num]/G[:,test_num])
		G[:, test_num] = G[:, test_num]*test_sdev
	G = G/np.std(G)
	return G

def standardize_ratio_between_expression_and_genotype(Y, G):
	num_tests = Y.shape[1]
	old_G = np.copy(G)
	for test_num in range(num_tests):
		rat = Y[:, test_num]/G[:,test_num]
		rat_sdev = np.std(rat)
		rat_mean = np.mean(rat)
		scaled_rat = (rat - rat_mean)/rat_sdev
		G[:, test_num] = Y[:, test_num]/scaled_rat
	G = G/np.std(G)
	return G

def standardize_columns(G):
	num_cols = G.shape[1]
	for col_num in range(num_cols):
		G[:, col_num] = (G[:, col_num] - np.mean(G[:, col_num]))/np.std(G[:, col_num])
	return G

def permute_donor(G, Z):
	G_perm = np.zeros(G.shape)
	num_samples = G.shape[0]
	num_tests = G.shape[1]
	num_donors = len(np.unique(Z))

	# First create donorXtest matrix
	G_donor = np.zeros((num_donors, num_tests))
	donor_num_to_sample_num = {}
	for donor_num in range(num_donors):
		donor_indices = np.where(Z==donor_num)[0]
		if len(donor_indices) < 1:
			print('Fatal assumption eroror')
			pdb.set_trace()
		sample_index = donor_indices[0]
		G_donor[donor_num, :] = G[sample_index, :]
		donor_num_to_sample_num[donor_num] = sample_index
	permy = np.random.permutation(range(num_donors))
	G_donor_perm = G_donor[permy,:]
	for sample_num in range(num_samples):
		G_perm[sample_num, :] = G_donor_perm[Z[sample_num], :]
	# generate sample level permy
	sample_permy = []
	for zz in Z:
		new_sample_num = donor_num_to_sample_num[permy[zz]]
		sample_permy.append(new_sample_num)
	sample_permy = np.asarray(sample_permy)
	if np.array_equal(G[sample_permy,:], G_perm) == False:
		print('assumption erororor')
	return G_perm, sample_permy

def get_resid_expression_from_vi_lmm(Y, G, Z, cov, lmm_root):
	F = np.load(lmm_root + 'temper_F.npy')
	alpha = np.load(lmm_root + 'temper_alpha.npy')
	C = np.load(lmm_root + 'temper_C.npy')

	N = Y.shape[0]
	T = Y.shape[1]

	# Random effects
	z_mapping = {}
	# Create mapping from grouping to index
	_, idx = np.unique(Z, return_index=True)
	unique_groups = np.asarray(Z)[np.sort(idx)]
	for i, label in enumerate(unique_groups):
		z_mapping[label] = i

	alpha_big_mu = np.zeros((N, T))
	for sample_num, z_label in enumerate(Z):
		alpha_big_mu[sample_num,:] = alpha[z_mapping[z_label], :]

	# Get residuals from VI model
	vi_pred_Y = np.dot(cov, C) + alpha_big_mu
	num_tests = G.shape[1]
	for test_num in range(num_tests):
		vi_pred_Y[:, test_num] = vi_pred_Y[:, test_num] + G[:, test_num]*F[test_num]

	Y_resid = Y - vi_pred_Y

	return Y_resid

def filter_RP_genes(test_names_df):
	num_genes = test_names_df.shape[0]
	valid_gene_indices = []
	for gene_num in range(num_genes):
		if test_names_df[gene_num,0].startswith('RP'):
			continue
		valid_gene_indices.append(gene_num)
	return np.asarray(valid_gene_indices)

def gene_correlated_with_used_expression_vectors(gene_vec, used_gene_vecs):
	gene_correlated = False
	for used_gene_vec in used_gene_vecs:
		corry = np.corrcoef(used_gene_vec, gene_vec)[0,1]
		if np.square(corry) > .1:
			gene_correlated = True
	return gene_correlated

def filter_RP_and_ind_genes(test_names_df, Y):
	num_genes = test_names_df.shape[0]
	valid_gene_indices = []
	used_expression_vectors = []

	for gene_num in np.random.permutation(num_genes):
		if test_names_df[gene_num,0].startswith('RP'):
			continue
		if gene_correlated_with_used_expression_vectors(Y[:, gene_num], used_expression_vectors):
			continue
		used_expression_vectors.append(Y[:, gene_num])
		valid_gene_indices.append(gene_num)
	return valid_gene_indices

def filter_RP_and_ind_genes_genos(test_names_df, Y, G):
	num_genes = test_names_df.shape[0]
	valid_gene_indices = []
	used_expression_vectors = []
	used_genotype_vectors = []

	counter = 0
	for gene_num in np.random.permutation(num_genes):
		counter = counter + 1
		if test_names_df[gene_num,0].startswith('RP'):
			continue
		if gene_correlated_with_used_expression_vectors(Y[:, gene_num], used_expression_vectors):
			continue
		if gene_correlated_with_used_expression_vectors(G[:, gene_num], used_genotype_vectors):
			continue
		used_expression_vectors.append(Y[:, gene_num])
		used_genotype_vectors.append(G[:, gene_num])
		valid_gene_indices.append(gene_num)
	return valid_gene_indices


def train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, covariate_file, num_latent_factors, output_root, model_name, lambda_v, variance_param, ard_variance_param, ratio_variance_standardization, permutation_type, warmup_iterations, round_genotype, data_filter, test_names_file, delta_elbo_threshold):
	# Load in expression data (dimension: num_samplesXnum_tests)
	Y = np.transpose(np.load(expression_training_file))
	# Load in Genotype data (dimension: num_samplesXnum_tests)
	G = np.transpose(np.load(genotype_training_file))
	# Load in sample overlap data
	Z,  num_individuals = load_in_sample_overlap_data(sample_overlap_file)
	# Load in covariates (dimension: num_samplesXnum_covariates)
	# We assume this covariate matrix has an intercept column
	cov = np.loadtxt(covariate_file)
	# Load in test names
	test_names_df = np.loadtxt(test_names_file, delimiter='\t',dtype=str)[1:,:]

	if round_genotype == "True":
		G = np.round(G)


	if data_filter == 'filter_RP':
		columns_filtered = filter_RP_genes(test_names_df)
		G = G[:, columns_filtered]
		Y = Y[:, columns_filtered]
		test_names_df = test_names_df[columns_filtered,:]
		np.savetxt(output_root + 'columns_filtered.txt', (columns_filtered), fmt="%s", delimiter='\t')
	elif data_filter == 'filter_RP_ind_genes':
		columns_filtered = filter_RP_and_ind_genes(test_names_df, Y)
		G = G[:, columns_filtered]
		Y = Y[:, columns_filtered]
		test_names_df = test_names_df[columns_filtered,:]
		np.savetxt(output_root + 'columns_filtered.txt', (columns_filtered), fmt="%s", delimiter='\t')
	elif data_filter == 'filter_RP_ind_genos':
		columns_filtered = filter_RP_and_ind_genes(test_names_df, G)
		G = G[:, columns_filtered]
		Y = Y[:, columns_filtered]
		test_names_df = test_names_df[columns_filtered,:]
		np.savetxt(output_root + 'columns_filtered.txt', (columns_filtered), fmt="%s", delimiter='\t')
	elif data_filter == 'filter_RP_ind_genes_genos':
		columns_filtered = filter_RP_and_ind_genes_genos(test_names_df, Y, G)
		G = G[:, columns_filtered]
		Y = Y[:, columns_filtered]
		test_names_df = test_names_df[columns_filtered,:]
		np.savetxt(output_root + 'columns_filtered.txt', (columns_filtered), fmt="%s", delimiter='\t')		





	if permutation_type == 'interaction_only':
		G_fe = np.copy(G)
		G, sample_permutation = permute_donor(G, Z)
		np.savetxt(output_root + 'sample_permutation.txt', (sample_permutation), fmt="%s", delimiter='\t')
	elif permutation_type == 'fixed_and_interaction':
		G, sample_permutation = permute_donor(G, Z)
		G_fe = np.copy(G)
		np.savetxt(output_root + 'sample_permutation.txt', (sample_permutation), fmt="%s", delimiter='\t')

	elif permutation_type == 'False':
		G_fe = np.copy(G)
	else:
		print('permutation type ' + permutation_type + ' currently not implemented')

	G_raw = np.copy(G)
	G = standardize_columns(G)
	G_fe = standardize_columns(G_fe)

	# Get number of samples, number of tests, number of individuals
	num_samples = Y.shape[0]
	num_tests = Y.shape[1]
	# Standardize variance of ratio between expression and genotype across tests
	if ratio_variance_standardization == 'True':
		G = standardize_variance_ratio_between_expression_and_genotype(Y, G)
		G_fe = standardize_variance_ratio_between_expression_and_genotype(Y, G_fe)



	#####################################
	# Run SURGE model
	#####################################
	if model_name == 'surge':
		re_boolean = True
		eqtl_vi = surge.surge_inference.SURGE_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, ard_alpha=ard_variance_param, ard_beta=ard_variance_param, max_iter=3000, gamma_v=lambda_v, warmup_iterations=warmup_iterations, re_boolean=re_boolean, delta_elbo_threshold=delta_elbo_threshold, verbose=True, output_root=output_root)
		eqtl_vi.fit(G=G, G_fe=G_fe, Y=Y, z=Z, cov=cov)
		
		# Save to output file
		np.savetxt(output_root + 'U_S.txt', (eqtl_vi.U_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'gamma_U.txt', eqtl_vi.gamma_U_alpha/eqtl_vi.gamma_U_beta, fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'V.txt', (eqtl_vi.V_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'alpha.txt', eqtl_vi.alpha_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'psi.txt', (eqtl_vi.psi_alpha/eqtl_vi.psi_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'elbo.txt', np.asarray(eqtl_vi.elbo), fmt="%s", delimiter='\n')
		np.savetxt(output_root + 'factor_genetic_pve.txt', (eqtl_vi.factor_genetic_pve), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'factor_pve.txt', (eqtl_vi.factor_pve), fmt="%s", delimiter='\t')





######################
# Command line args
######################
expression_training_file = sys.argv[1]
genotype_training_file = sys.argv[2]
covariate_file = sys.argv[3]
sample_overlap_file = sys.argv[4]
num_latent_factors = int(sys.argv[5])
lambda_v = float(sys.argv[6])
model_name = sys.argv[7]
seed = int(sys.argv[8])
output_root = sys.argv[9]
variance_param = float(sys.argv[10])
ard_variance_param = float(sys.argv[11])
ratio_variance_standardization = sys.argv[12]
permutation_type = sys.argv[13]
warmup_iterations = int(sys.argv[14])
round_genotype = sys.argv[15]
data_filter = sys.argv[16]
test_names_file = sys.argv[17]
delta_elbo_threshold = float(sys.argv[18])


np.random.seed(seed)


train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, covariate_file, num_latent_factors, output_root, model_name, lambda_v, variance_param, ard_variance_param, ratio_variance_standardization, permutation_type, warmup_iterations, round_genotype, data_filter, test_names_file, delta_elbo_threshold)


