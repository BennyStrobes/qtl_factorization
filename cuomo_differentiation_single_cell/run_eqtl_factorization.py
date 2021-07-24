import numpy as np 
import os
import sys
import pdb
import eqtl_factorization_vi_ard



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
	for donor_num in range(num_donors):
		donor_indices = np.where(Z==donor_num)[0]
		if len(donor_indices) < 1:
			print('Fatal assumption eroror')
			pdb.set_trace()
		sample_index = donor_indices[0]
		G_donor[donor_num, :] = G[sample_index, :]
	permy = np.random.permutation(range(num_donors))
	G_donor_perm = G_donor[permy,:]
	for sample_num in range(num_samples):
		G_perm[sample_num, :] = G_donor_perm[Z[sample_num], :]
	return G_perm

def train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, covariate_file, num_latent_factors, output_root, model_name, lambda_v, variance_param, ard_variance_param, ratio_variance_standardization, permutation_type, warmup_iterations):
	# Load in expression data (dimension: num_samplesXnum_tests)
	Y = np.transpose(np.load(expression_training_file))
	# Load in Genotype data (dimension: num_samplesXnum_tests)
	G = np.transpose(np.load(genotype_training_file))
	# Load in sample overlap data
	Z,  num_individuals = load_in_sample_overlap_data(sample_overlap_file)
	# Load in covariates (dimension: num_samplesXnum_covariates)
	# We assume this covariate matrix has an intercept column
	cov = np.loadtxt(covariate_file)

	if permutation_type == 'interaction_only':
		G_fe = np.copy(G)
		G = permute_donor(G, Z)
	elif permutation_type == 'fixed_and_interaction':
		G = permute_donor(G, Z)
		G_fe = np.copy(G)
	elif permutation_type == 'False':
		G_fe = np.copy(G)
	else:
		print('permutation type ' + permutation_type + ' currently not implemented')

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
	if model_name == 'eqtl_factorization_vi_ard':
		eqtl_vi = eqtl_factorization_vi_ard.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, ard_alpha=ard_variance_param, ard_beta=ard_variance_param, max_iter=400, gamma_v=lambda_v, warmup_iterations=warmup_iterations, output_root=output_root)
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


np.random.seed(seed)

train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, covariate_file, num_latent_factors, output_root, model_name, lambda_v, variance_param, ard_variance_param, ratio_variance_standardization, permutation_type, warmup_iterations)


