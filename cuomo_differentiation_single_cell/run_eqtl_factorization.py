import numpy as np 
import os
import sys
import pdb
import eqtl_factorization_vi_ard
import eqtl_factorization_vi_ard_heteroskedastic
import eqtl_factorization_vi_ard_no_F
import eqtl_factorization_vi_ard_alt_init
import eqtl_factorization_vi_no_factorization
from pymer4.models import Lmer
import pandas as pd
import eqtl_factorization_no_factorization
import eqtl_factorization_vi_no_factorization_lmm_init
import eqtl_factorization_vi_fixed_per_test
import eqtl_factorization_vi_ard_test_specific_tau_prior
import eqtl_factorization_vi_ard_fixed_test_specific_tau
import eqtl_factorization_vi_ard_environmental_fixed_effect
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

def extract_residuals_from_standard_eqtl_model(Y, G, cov, z):
	num_tests = Y.shape[1]
	F_betas = []
	C_betas = []
	residuals = []
	N = Y.shape[0]
	model_eq = 'y ~ g'
	for cov_num in range(cov.shape[1]):
		model_eq = model_eq + ' + x' + str(cov_num)
	model_eq = model_eq + ' + (1|z)'
	# 119, 103
	#for test_number in range(num_tests):
	vec1 = []
	vec2 = []
	for test_number in [0,1,2,3,4,5,6,7,8,9,10,103]:
		print(test_number)
		y_vec = Y[:,test_number]
		g_vec = G[:,test_number]
		np.savetxt('y.txt', y_vec, fmt="%s", delimiter='\t')
		np.savetxt('g.txt', g_vec, fmt="%s", delimiter='\t')
		np.savetxt('cov.txt', cov, fmt="%s", delimiter='\t')
		np.savetxt('z.txt', z, fmt="%s", delimiter='\t')
		os.system('Rscript lmm.R')
		#dd = {'y':y_vec, 'z':z, 'g':g_vec}
		#num_covs = cov.shape[1]
		#for cov_num in range(num_covs):
			#dd['x' + str(cov_num)] = cov[:, cov_num]
		#df = pd.DataFrame(dd)
		#model = Lmer(model_eq, data=df)
		#model.fit()
		#residuals.append(model.residuals)
		#resid = model.residuals
		#ss_ratio = np.sum(np.square(resid/g_vec))
		#ss2 = np.sum(np.square(resid))*ef_tau/2.0
		model_residuals = np.loadtxt('lmm_residuals.txt')
		residual_variance = np.loadtxt('residual_variance.txt')*1.0
		ef_tau = 1.0/residual_variance
		print(residual_variance)

		ss2 = np.sum(np.square(model_residuals))*ef_tau/2.0 + ((N/2.0)*np.log(2.0*np.pi*(1.0/ef_tau)))
		ss_ratio = np.log(np.sum(np.square(model_residuals/g_vec)))
		vec1.append(ss2)
		vec2.append(ss_ratio)
		#print(np.mean(model.residuals/g_vec)/np.std(model.residuals/g_vec))
		#print('\n')
	pdb.set_trace()
	residuals = np.transpose(np.asarray(residuals))
	return residuals

def remove_outlier_tests(Y, G, G_fe):
	means = np.mean(Y/G,axis=0)
	upper_bound = np.mean(means) + 2.0*np.std(means)
	lower_bound = np.mean(means) - 2.0*np.std(means)
	new_test_indices = (means > lower_bound)*(means < upper_bound)
	return Y[:, new_test_indices], G[:, new_test_indices], G_fe[:, new_test_indices]

def get_lmm_residual_expression(Y, G, G_fe, Z, model_covariates):
	# LOAD IN VI LMM model
	file_stem = '/work-zfs/abattle4/bstrober/qtl_factorization/cuomo_differentiation_single_cell/eqtl_factorization_results/'
	#file_stem = file_stem + 'eqtl_factorization_standard_eqtl_scran_1000_hvg_eqtl_factorization_vi_no_factorization_lmm_init_results_k_init_1_seed_1_warmup_1000_ratio_variance_std_True_permute_False_no_F_prior_temper_'
	file_stem = file_stem + 'eqtl_factorization_standard_eqtl_scran_1000_hvg_eqtl_factorization_vi_no_factorization_results_k_init_1_seed_1_warmup_1000_ratio_variance_std_False_permute_False_no_F_prior_temper_'

	eqtl_factorization_C_file = file_stem + 'C.npy'
	vi_C = np.load(eqtl_factorization_C_file)

	eqtl_factorization_F_file = file_stem + 'F.npy'
	vi_F = np.load(eqtl_factorization_F_file)

	eqtl_factorization_alpha_file = file_stem + 'alpha.npy'
	vi_alpha = np.load(eqtl_factorization_alpha_file)

	eqtl_factorization_tau_file = file_stem + 'tau.npy'
	vi_tau = np.load(eqtl_factorization_tau_file)

	eqtl_factorization_psi_file = file_stem + 'psi.npy'
	vi_psi = np.load(eqtl_factorization_psi_file)

	# Get donor means from each sample for vi model
	num_samples = len(Z)
	vi_alpha_big = []
	for sample_num in range(num_samples):
		individual_corresponding_to_sample = Z[sample_num]
		vi_alpha_big.append(vi_alpha[individual_corresponding_to_sample,:])
	vi_alpha_big = np.asarray(vi_alpha_big)

	G_temp = np.copy(G)

	# Get residuals from VI model
	vi_pred_Y = np.dot(model_covariates, vi_C) + vi_alpha_big
	num_tests = G_temp.shape[1]
	for test_num in range(num_tests):
		vi_pred_Y[:, test_num] = vi_pred_Y[:, test_num] + G_temp[:, test_num]*vi_F[test_num]  
	vi_resid_Y = Y - vi_pred_Y

	return vi_resid_Y

def remove_ratio_outliers(Y, G, G_fe, Z, model_covariates):
	# LOAD IN VI LMM model
	file_stem = '/work-zfs/abattle4/bstrober/qtl_factorization/cuomo_differentiation_single_cell/eqtl_factorization_results/'
	#file_stem = file_stem + 'eqtl_factorization_standard_eqtl_scran_1000_hvg_eqtl_factorization_vi_no_factorization_lmm_init_results_k_init_1_seed_1_warmup_1000_ratio_variance_std_True_permute_False_no_F_prior_temper_'
	file_stem = file_stem + 'eqtl_factorization_standard_eqtl_scran_1000_hvg_eqtl_factorization_vi_no_factorization_results_k_init_1_seed_1_warmup_1000_ratio_variance_std_False_permute_False_no_F_prior_temper_'

	eqtl_factorization_C_file = file_stem + 'C.npy'
	vi_C = np.load(eqtl_factorization_C_file)

	eqtl_factorization_F_file = file_stem + 'F.npy'
	vi_F = np.load(eqtl_factorization_F_file)

	eqtl_factorization_alpha_file = file_stem + 'alpha.npy'
	vi_alpha = np.load(eqtl_factorization_alpha_file)

	eqtl_factorization_tau_file = file_stem + 'tau.npy'
	vi_tau = np.load(eqtl_factorization_tau_file)

	eqtl_factorization_psi_file = file_stem + 'psi.npy'
	vi_psi = np.load(eqtl_factorization_psi_file)

	# Get donor means from each sample for vi model
	num_samples = len(Z)
	vi_alpha_big = []
	for sample_num in range(num_samples):
		individual_corresponding_to_sample = Z[sample_num]
		vi_alpha_big.append(vi_alpha[individual_corresponding_to_sample,:])
	vi_alpha_big = np.asarray(vi_alpha_big)

	G_temp = np.copy(G)

	# Get residuals from VI model
	vi_pred_Y = np.dot(model_covariates, vi_C) + vi_alpha_big
	num_tests = G_temp.shape[1]
	for test_num in range(num_tests):
		vi_pred_Y[:, test_num] = vi_pred_Y[:, test_num] + G_temp[:, test_num]*vi_F[test_num]  
	vi_resid_Y = Y - vi_pred_Y
	ratio = (np.var(vi_resid_Y/G_temp,axis=0)/np.var(Y/G_temp,axis=0))/np.var(vi_resid_Y,axis=0)

	abs_diff = np.abs(ratio-1.0)

	new_test_indices = abs_diff < .01

	return Y[:, new_test_indices], G[:, new_test_indices], G_fe[:, new_test_indices], vi_resid_Y[:, new_test_indices]

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

	G_raw = np.copy(G)

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

	remove_ratio_outliers_bool = 'False'
	if remove_ratio_outliers_bool == 'True':
		Y, G, G_fe, Y_resid = remove_ratio_outliers(Y, G, G_fe, Z, cov)
	elif remove_ratio_outliers_bool == 'False':
		Y_resid = get_lmm_residual_expression(Y, G, G_fe, Z, cov)

	# Get number of samples, number of tests, number of individuals
	num_samples = Y.shape[0]
	num_tests = Y.shape[1]
	# Standardize variance of ratio between expression and genotype across tests
	if ratio_variance_standardization == 'True':
		G = standardize_variance_ratio_between_expression_and_genotype(Y, G)
		G_fe = standardize_variance_ratio_between_expression_and_genotype(Y, G_fe)
	elif ratio_variance_standardization == 'True_residuals':
		G = standardize_variance_ratio_between_expression_and_genotype(Y_resid, G)
		G_fe = standardize_variance_ratio_between_expression_and_genotype(Y_resid, G_fe)	

	'''
	valid_indices = np.asarray([True]*Y.shape[1])
	corrz1_arr = []
	corrz2_arr = []
	for test_num in range(num_tests):
		corrz1 = np.abs(np.corrcoef(np.abs(Y[:, test_num]), G[:, test_num])[0,1])
		corrz2 = np.abs(np.corrcoef(np.abs(Y_resid[:, test_num]), G[:, test_num])[0,1])
		corrz1_arr.append(corrz1)
		corrz2_arr.append(corrz2)
		if corrz1 > .2:
			valid_indices[test_num] = False
	Y = Y[:, valid_indices]
	G = G[:, valid_indices]
	G_fe = G_fe[:, valid_indices]
	'''


	#####################################
	# Run SURGE model
	#####################################
	if model_name == 'eqtl_factorization_vi_ard':
		eqtl_vi = eqtl_factorization_vi_ard.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, ard_alpha=ard_variance_param, ard_beta=ard_variance_param, max_iter=600, gamma_v=lambda_v, warmup_iterations=warmup_iterations, output_root=output_root)
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
	elif model_name == 'eqtl_factorization_vi_ard_environmental_fixed_effect':
		eqtl_vi = eqtl_factorization_vi_ard_environmental_fixed_effect.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, ard_alpha=ard_variance_param, ard_beta=ard_variance_param, max_iter=600, gamma_v=lambda_v, warmup_iterations=warmup_iterations, output_root=output_root)
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
	elif model_name == 'eqtl_factorization_vi_ard_heteroskedastic':
		eqtl_vi = eqtl_factorization_vi_ard_heteroskedastic.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, ard_alpha=ard_variance_param, ard_beta=ard_variance_param, max_iter=600, gamma_v=lambda_v, warmup_iterations=warmup_iterations, output_root=output_root)
		eqtl_vi.fit(G=G, G_fe=G_fe, G_raw=G_raw, Y=Y, z=Z, cov=cov)

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
	elif model_name == 'eqtl_factorization_vi_ard_test_specific_tau_prior':
		tau_prior_file = '/work-zfs/abattle4/bstrober/qtl_factorization/cuomo_differentiation_single_cell/eqtl_factorization_results/eqtl_factorization_standard_eqtl_scran_1000_hvg_eqtl_factorization_vi_no_factorization_results_k_init_1_seed_1_warmup_5000_ratio_variance_std_True_permute_False_temper_tau.npy'
		tau_prior = np.load(tau_prior_file)
		eqtl_vi = eqtl_factorization_vi_ard_test_specific_tau_prior.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, tau_prior=tau_prior, ard_alpha=ard_variance_param, ard_beta=ard_variance_param, max_iter=600, gamma_v=lambda_v, warmup_iterations=warmup_iterations, output_root=output_root)
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
	elif model_name == 'eqtl_factorization_vi_ard_fixed_test_specific_tau':
		tau_prior_file = '/work-zfs/abattle4/bstrober/qtl_factorization/cuomo_differentiation_single_cell/eqtl_factorization_results/eqtl_factorization_standard_eqtl_scran_1000_hvg_eqtl_factorization_vi_no_factorization_results_k_init_1_seed_1_warmup_5000_ratio_variance_std_True_permute_False_temper_tau.npy'
		tau_prior = np.load(tau_prior_file)
		eqtl_vi = eqtl_factorization_vi_ard_fixed_test_specific_tau.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, tau_prior=tau_prior, ard_alpha=ard_variance_param, ard_beta=ard_variance_param, max_iter=600, gamma_v=lambda_v, warmup_iterations=warmup_iterations, output_root=output_root)
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
	elif model_name == 'eqtl_factorization_vi_ard_no_F':
		eqtl_vi = eqtl_factorization_vi_ard_no_F.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, ard_alpha=ard_variance_param, ard_beta=ard_variance_param, max_iter=600, gamma_v=lambda_v, warmup_iterations=warmup_iterations, output_root=output_root)
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
	elif model_name == 'eqtl_factorization_vi_ard_alt_init':
		eqtl_vi = eqtl_factorization_vi_ard_alt_init.EQTL_FACTORIZATION_VI(K=num_latent_factors, alpha=variance_param, beta=variance_param, ard_alpha=ard_variance_param, ard_beta=ard_variance_param, max_iter=400, gamma_v=lambda_v, warmup_iterations=warmup_iterations, output_root=output_root)
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
	elif model_name == 'eqtl_factorization_vi_no_factorization':
		eqtl_vi = eqtl_factorization_vi_no_factorization.EQTL_FACTORIZATION_VI(alpha=variance_param, beta=variance_param, max_iter=200, output_root=output_root)
		eqtl_vi.fit(G=G, G_fe=G_fe, Y=Y, z=Z, cov=cov)
		# Save to output file
		np.savetxt(output_root + 'F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'alpha.txt', eqtl_vi.alpha_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'psi.txt', (eqtl_vi.psi_alpha/eqtl_vi.psi_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'elbo.txt', np.asarray(eqtl_vi.elbo), fmt="%s", delimiter='\n')
	elif model_name == 'eqtl_factorization_no_factorization':
		eqtl_vi = eqtl_factorization_no_factorization.EQTL_FACTORIZATION(output_root=output_root)
		eqtl_vi.fit(G=G, G_fe=G_fe, Y=Y, z=Z, cov=cov)
		# Save to output file
		np.savetxt(output_root + 'F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'alpha.txt', eqtl_vi.alpha_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'tau.txt', (eqtl_vi.tau), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'psi.txt', (eqtl_vi.psi), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
	elif model_name == 'eqtl_factorization_vi_no_factorization_lmm_init':
		eqtl_vi = eqtl_factorization_vi_no_factorization_lmm_init.EQTL_FACTORIZATION_VI(alpha=variance_param, beta=variance_param, max_iter=600, output_root=output_root)
		eqtl_vi.fit(G=G, G_fe=G_fe, Y=Y, z=Z, cov=cov)
		# Save to output file
		np.savetxt(output_root + 'F.txt', (eqtl_vi.F_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'alpha.txt', eqtl_vi.alpha_mu, fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'tau.txt', (eqtl_vi.tau_alpha/eqtl_vi.tau_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'psi.txt', (eqtl_vi.psi_alpha/eqtl_vi.psi_beta), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'C.txt', (eqtl_vi.C_mu), fmt="%s", delimiter='\t')
		np.savetxt(output_root + 'elbo.txt', np.asarray(eqtl_vi.elbo), fmt="%s", delimiter='\n')
	elif model_name == 'eqtl_factorization_vi_fixed_per_test':
		delta_elbo_threshold = 1e-100
		max_iter = 40
		elbos = []
		output_file = output_root + 'per_test_elbos.txt'
		t = open(output_file,'w')
		#for test_num in range(Y.shape[1]):
		for test_num in [103,119]:
			eqtl_vi = eqtl_factorization_vi_fixed_per_test.EQTL_FACTORIZATION_VI(test_num=test_num, alpha=variance_param, beta=variance_param, max_iter=max_iter, delta_elbo_threshold=delta_elbo_threshold, output_root=output_root)
			eqtl_vi.fit(G=G, G_fe=G_fe, Y=Y, z=Z, cov=cov)
			t.write('\t'.join(np.asarray(eqtl_vi.elbo).astype(str)) + '\n')
			t.flush()
		t.close()
	elif model_name == 'eqtl_factorization_vi_null_fixed_per_test':
		delta_elbo_threshold = 1e-100
		max_iter = 40
		test_num = 0
		output_file = output_root + 'per_test_elbos.txt'
		t = open(output_file,'w')
		for test_num in range(Y.shape[1]):
			eqtl_vi = eqtl_factorization_vi_fixed_per_test.EQTL_FACTORIZATION_VI(test_num=test_num, alpha=variance_param, beta=variance_param, max_iter=max_iter, delta_elbo_threshold=delta_elbo_threshold, output_root=output_root)
			eqtl_vi.fit(G=np.ones(G.shape), G_fe=G_fe, Y=Y, z=Z, cov=cov)
			t.write('\t'.join(np.asarray(eqtl_vi.elbo).astype(str)) + '\n')
			t.flush()
		t.close()








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


