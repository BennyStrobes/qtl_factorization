import numpy as np 
import os
import sys
import pdb
import eqtl_factorization_vi
import eqtl_factorization_vi_ard
import eqtl_factorization_pca
import eqtl_factorization_vi_bernoulli_loading
import eqtl_factorization_vi_fixed_residual_variance
import eqtl_factorization_vi_lda
import eqtl_factorization_vi_lda_gaussian_factors
import eqtl_factorization_vi_hdp
import eqtl_factorization_vi_ard_factors_gaussian_loadings
import eqtl_factorization_vi_gaussian_factors_gaussian_loadings
import eqtl_factorization_vi_ard_permute_k
import eqtl_factorization_vi_ard_no_re



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

def train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, covariate_file, num_latent_factors, output_root, model_name, lambda_v, variance_param, ratio_variance_standardization, permutation_type):
	# Load in expression data (dimension: num_samplesXnum_tests)
	Y = np.transpose(np.load(expression_training_file))
	# Load in Genotype data (dimension: num_samplesXnum_tests)
	G = np.transpose(np.load(genotype_training_file))
	# Load in sample overlap data
	Z,  num_individuals = load_in_sample_overlap_data(sample_overlap_file)
	# Load in covariates (dimension: num_samplesXnum_covariates)
	# We assume this covariate matrix has an intercept column
	cov = np.loadtxt(covariate_file)

	# Get number of samples, number of tests, number of individuals
	num_samples = Y.shape[0]
	num_tests = Y.shape[1]
	# Standardize variance of ratio between expression and genotype across tests
	if ratio_variance_standardization == 'True':
		G = standardize_variance_ratio_between_expression_and_genotype(Y, G)
	elif ratio_variance_standardization == 'Alt':
		G = G*100.0

	#####################################
	# Run model
	#####################################
	if model_name == 'eqtl_factorization_vi':
		eqtl_vi = eqtl_factorization_vi.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, a=1, b=1, max_iter=800, gamma_v=lambda_v, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, z=Z, cov=cov)
		# Order and Filter Factors
		theta_U = eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a)
		ordered_indices = np.argsort(-theta_U)
		num_indices = sum(theta_U > .01)
		ordered_filtered_indices = ordered_indices[:num_indices]
		# Save to output file
		np.savetxt(output_root + 'V.txt', (eqtl_vi.V_mu)[ordered_filtered_indices, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'theta_U.txt', theta_U[ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'S.txt', (eqtl_vi.S_U)[:, ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_ard':
		eqtl_vi = eqtl_factorization_vi_ard.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, a=1, b=1, max_iter=400, gamma_v=lambda_v, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, z=Z, cov=cov)
		# Order and Filter Factors
		theta_U = eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a)
		ordered_indices = np.argsort(-theta_U)
		num_indices = sum(theta_U > .01)
		ordered_filtered_indices = ordered_indices[:num_indices]
		# Save to output file
		np.savetxt(output_root + 'V.txt', (eqtl_vi.V_mu)[ordered_filtered_indices, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'theta_U.txt', theta_U[ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'S.txt', (eqtl_vi.S_U)[:, ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_ard_no_re':
		eqtl_vi = eqtl_factorization_vi_ard_no_re.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, a=1, b=1, max_iter=400, gamma_v=lambda_v, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, z=Z, cov=cov)
		# Order and Filter Factors
		theta_U = eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a)
		ordered_indices = np.argsort(-theta_U)
		num_indices = sum(theta_U > .01)
		ordered_filtered_indices = ordered_indices[:num_indices]
		# Save to output file
		np.savetxt(output_root + 'V.txt', (eqtl_vi.V_mu)[ordered_filtered_indices, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'theta_U.txt', theta_U[ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'S.txt', (eqtl_vi.S_U)[:, ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_ard_permute_k':
		eqtl_vi = eqtl_factorization_vi_ard_permute_k.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, a=1, b=1, max_iter=200, gamma_v=lambda_v, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, z=Z, cov=cov)
		# Order and Filter Factors
		theta_U = eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a)
		ordered_indices = np.argsort(-theta_U)
		num_indices = sum(theta_U > .01)
		ordered_filtered_indices = ordered_indices[:num_indices]
		# Save to output file
		np.savetxt(output_root + 'V.txt', (eqtl_vi.V_mu)[ordered_filtered_indices, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'theta_U.txt', theta_U[ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'S.txt', (eqtl_vi.S_U)[:, ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_ard_factors_gaussian_loadings':
		eqtl_vi = eqtl_factorization_vi_ard_factors_gaussian_loadings.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, a=1, b=1, max_iter=800, gamma_v=lambda_v, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, z=Z, cov=cov)
		# Order and Filter Factors
		theta_U = eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a)
		ordered_indices = np.argsort(-theta_U)
		num_indices = sum(theta_U > .01)
		ordered_filtered_indices = ordered_indices[:num_indices]
		# Save to output file
		np.savetxt(output_root + 'V.txt', (eqtl_vi.V_mu)[ordered_filtered_indices, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'theta_U.txt', theta_U[ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'S.txt', (eqtl_vi.S_U)[:, ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_ard_factors_gaussian_loadings_variance_weighted':
		eqtl_vi = eqtl_factorization_vi_ard_factors_gaussian_loadings.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, a=1, b=1, max_iter=800, gamma_v=lambda_v, output_root=output_root)
		num_tests = Y.shape[1]
		for test_num in range(num_tests):
			test_sdev = np.std(Y[:, test_num]/G[:,test_num])
			G[:, test_num] = G[:, test_num]*test_sdev
		G = G/np.std(G)
		eqtl_vi.fit(G=G, Y=Y, z=Z, cov=cov)
		# Order and Filter Factors
		theta_U = eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a)
		ordered_indices = np.argsort(-theta_U)
		num_indices = sum(theta_U > .01)
		ordered_filtered_indices = ordered_indices[:num_indices]
		# Save to output file
		np.savetxt(output_root + 'V.txt', (eqtl_vi.V_mu)[ordered_filtered_indices, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'theta_U.txt', theta_U[ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'S.txt', (eqtl_vi.S_U)[:, ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_gaussian_factors_gaussian_loadings':
		eqtl_vi = eqtl_factorization_vi_gaussian_factors_gaussian_loadings.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, a=1, b=1, max_iter=800, gamma_v=lambda_v, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, z=Z, cov=cov)
		# Order and Filter Factors
		theta_U = eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a)
		ordered_indices = np.argsort(-theta_U)
		num_indices = sum(theta_U > .01)
		ordered_filtered_indices = ordered_indices[:num_indices]
		# Save to output file
		np.savetxt(output_root + 'V.txt', (eqtl_vi.V_mu)[ordered_filtered_indices, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'theta_U.txt', theta_U[ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'S.txt', (eqtl_vi.S_U)[:, ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_lda':
		eqtl_vi = eqtl_factorization_vi_lda.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, a=1, b=1, max_iter=800, delta_0=.5, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, z=Z, cov=cov)
		# Order and Filter Factors
		theta_U = eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a)
		ordered_indices = np.argsort(-theta_U)
		num_indices = sum(theta_U > .01)
		ordered_filtered_indices = ordered_indices[:num_indices]
		# Save to output file
		np.savetxt(output_root + 'V.txt', (eqtl_vi.V_mu)[ordered_filtered_indices, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'theta_U.txt', theta_U[ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'S.txt', (eqtl_vi.S_U)[:, ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_lda_gaussian_factors':
		eqtl_vi = eqtl_factorization_vi_lda_gaussian_factors.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, a=1, b=1, max_iter=800, delta_0=.5, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, z=Z, cov=cov)
		# Order and Filter Factors
		theta_U = eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a)
		ordered_indices = np.argsort(-theta_U)
		num_indices = sum(theta_U > .01)
		ordered_filtered_indices = ordered_indices[:num_indices]
		# Save to output file
		np.savetxt(output_root + 'V.txt', (eqtl_vi.V_mu)[ordered_filtered_indices, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'theta_U.txt', theta_U[ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'S.txt', (eqtl_vi.S_U)[:, ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_hdp':
		eqtl_vi = eqtl_factorization_vi_hdp.EQTL_FACTORIZATION_VI(K=num_latent_factors, J=num_latent_factors-5, alpha=variance_param, beta=variance_param, a=1, b=1, max_iter=800, delta_0=1.0, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, z=Z, cov=cov)
		# Order and Filter Factors
		theta_U = eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a)
		ordered_indices = np.argsort(-theta_U)
		num_indices = sum(theta_U > .01)
		ordered_filtered_indices = ordered_indices[:num_indices]
		# Save to output file
		np.savetxt(output_root + 'V.txt', (eqtl_vi.V_mu)[ordered_filtered_indices, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'theta_U.txt', theta_U[ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'S.txt', (eqtl_vi.S_U)[:, ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	if model_name == 'eqtl_factorization_vi_fixed_residual_variance':
		eqtl_vi = eqtl_factorization_vi_fixed_residual_variance.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, a=1, b=1, max_iter=800, gamma_v=lambda_v, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, z=Z, cov=cov)
		# Order and Filter Factors
		theta_U = eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a)
		ordered_indices = np.argsort(-theta_U)
		num_indices = sum(theta_U > .01)
		ordered_filtered_indices = ordered_indices[:num_indices]
		# Save to output file
		np.savetxt(output_root + 'V.txt', (eqtl_vi.V_mu)[ordered_filtered_indices, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'theta_U.txt', theta_U[ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'tau.txt', [eqtl_vi.tau_alpha/eqtl_vi.tau_beta], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'S.txt', (eqtl_vi.S_U)[:, ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	elif model_name == 'eqtl_factorization_pca':
		eqtl_pca = eqtl_factorization_pca.EQTL_FACTORIZATION_PCA(K=num_latent_factors, output_root=output_root)
		eqtl_pca.fit(G=G, Y=Y, z=Z, cov=cov)
	elif model_name == 'eqtl_factorization_vi_bernoulli_loading':
		eqtl_vi = eqtl_factorization_vi_bernoulli_loading.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, a=1, b=1, max_iter=300, gamma_v=lambda_v, output_root=output_root)
		eqtl_vi.fit(G=G, Y=Y, z=Z, cov=cov)
		# Order and Filter Factors
		theta_U = eqtl_vi.theta_U_a/(eqtl_vi.theta_U_b + eqtl_vi.theta_U_a)
		ordered_indices = np.argsort(-theta_U)
		num_indices = sum(theta_U > .01)
		ordered_filtered_indices = ordered_indices[:num_indices]
		# Save to output file
		np.savetxt(output_root + 'V.txt', (eqtl_vi.V_mu)[ordered_filtered_indices, :], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'U_S.txt', (eqtl_vi.U_mu*eqtl_vi.S_U)[:,ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'theta_U.txt', theta_U[ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'S.txt', (eqtl_vi.S_U)[:, ordered_filtered_indices], fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')










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
ratio_variance_standardization = sys.argv[11]
permutation_type = sys.argv[12]


np.random.seed(seed)


train_eqtl_factorization_model(sample_overlap_file, expression_training_file, genotype_training_file, covariate_file, num_latent_factors, output_root, model_name, lambda_v, variance_param, ratio_variance_standardization, permutation_type)


