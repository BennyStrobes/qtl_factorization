import numpy as np 
import os
import sys
import pdb
import eqtl_factorization_ard_no_alpha
from sklearn.linear_model import LinearRegression
from statsmodels.stats import diagnostic

def get_genotype_vec(num_samples, af):
	allele1 = np.random.binomial(n=1, p=af,size=num_samples)
	allele2 = np.random.binomial(n=1, p=af,size=num_samples)
	genotype_vec = allele1 + allele2
	return genotype_vec

def generate_null_eqtl_factorization_data_with_no_sample_repeat_and_first_test_has_heteroskedasticity(num_samples, num_tests, test_0_af, version):
	# Simulate genotype
	G = np.zeros((num_samples, num_tests))
	G_raw = np.zeros((num_samples, num_tests))
	mafs = []
	for test_num in range(num_tests):
		if test_num == 0:
			af = np.copy(test_0_af)
		else:
			af = np.random.uniform(.05,.95)
		# Get genotype for this particular snp
		genotype_vec = get_genotype_vec(num_samples, af)
		if np.mean(genotype_vec) == 1.0:
			genotype_vec = get_genotype_vec(num_samples, af)
			if np.mean(genotype_vec) == 1.0:
				genotype_vec = get_genotype_vec(num_samples, af)
				if np.mean(genotype_vec) == 1.0:
					genotype_vec = get_genotype_vec(num_samples, af)
		af = np.sum(genotype_vec)/(2.0*len(genotype_vec))
		if af > .5:
			maf = 1.0 - af
		else:
			maf = af
		mafs.append(maf)
		G_raw[:, test_num] = genotype_vec
		G[:, test_num] = (genotype_vec - np.mean(genotype_vec))/np.std(genotype_vec)
	mafs = np.asarray(mafs)
	predicted_mean = np.zeros((num_samples, num_tests))
	betas = []
	Y = np.zeros((num_samples, num_tests))
	for test_num in range(num_tests):
		if test_num == 0 and version == 'heteroskedasticity':
			beta = np.random.normal(loc=0.0, scale=.1)
			predicted_mean[:, test_num] = predicted_mean[:, test_num] + G[:, test_num]*beta
			indices0 = G_raw[:,test_num] == 0.0
			indices1 = G_raw[:,test_num] == 1.0
			indices2 = G_raw[:,test_num] == 2.0

			Y[indices0,test_num] = np.random.normal(predicted_mean[indices0,test_num], scale=.1)
			Y[indices1,test_num] = np.random.normal(predicted_mean[indices1,test_num], scale=1)
			Y[indices2,test_num] = np.random.normal(predicted_mean[indices2,test_num], scale=5)
		else:
			beta = np.random.normal(loc=0.0, scale=.1)
			predicted_mean[:, test_num] = predicted_mean[:, test_num] + G[:, test_num]*beta
			Y[:, test_num] = np.random.normal(predicted_mean[:,test_num])
		betas.append(beta)
	for test_num in range(num_tests):
		Y[:, test_num] = (Y[:, test_num] - np.mean(Y[:, test_num]))/np.std(Y[:, test_num])
	return Y, G, np.asarray(betas), mafs, G_raw


def get_residuals(Y, G):
	num_tests = Y.shape[1]
	num_samples = Y.shape[0]
	learned_betas = []
	resids = []
	corrz = []
	for test_number in range(num_tests):
		y_vec = Y[:,test_number]
		g_vec = G[:,test_number]
		X = np.hstack((np.transpose(np.asmatrix(g_vec)), np.ones((num_samples, 1))))
		reg = LinearRegression(fit_intercept=False).fit(X, np.transpose(np.asmatrix(y_vec)))
		pred_y = reg.predict(X)[:,0]
		resid_y = y_vec - pred_y
		learned_beta = reg.coef_[0][0]
		learned_betas.append(learned_beta)
		resids.append(resid_y)
		corry = np.corrcoef(np.abs(resid_y), g_vec)[0,1]
		corrz.append(corry)
	return np.transpose(np.asarray(resids)), np.asarray(learned_betas), np.asarray(corrz)

####################
# Command line args
####################
num_samples = int(sys.argv[1])
num_tests = int(sys.argv[2])
af = float(sys.argv[3])
seed = int(sys.argv[4])
version = sys.argv[5]
output_file = sys.argv[6]

np.random.seed(seed)



Y, G, betas, mafs, G_raw = generate_null_eqtl_factorization_data_with_no_sample_repeat_and_first_test_has_heteroskedasticity(num_samples, num_tests, af, version)
Y_resid, beta_estimated, corrz = get_residuals(Y, G)
Y_pred = Y-Y_resid

heteroskedastic_pvalz = []

for test_iter in range(Y_pred.shape[1]):
	exog = np.hstack((np.ones((len(Y_pred[:,test_iter]),1)), np.transpose(np.asmatrix(Y_pred[:,test_iter]))))
	_, pval, __, f_pval = diagnostic.het_white(Y_resid[:,test_iter], exog)
	heteroskedastic_pvalz.append(f_pval)
heteroskedastic_pvalz = np.asarray(heteroskedastic_pvalz)

# Dummy variable
Z = np.ones(num_samples)
cov = np.ones((num_samples, 1))



eqtl_vi = eqtl_factorization_ard_no_alpha.EQTL_FACTORIZATION_VI(K=3, alpha=1e-3, beta=1e-3, a=1, b=1, max_iter=20, gamma_v=1)
eqtl_vi.fit(G=G, Y=Y, z=Z, cov=cov)

test_tau = (eqtl_vi.tau_alpha/eqtl_vi.tau_beta)[0]

strongest_component = np.argmin(eqtl_vi.gamma_U_alpha/eqtl_vi.gamma_U_beta)

if np.argmax(np.abs(eqtl_vi.V_mu[strongest_component,:])) == 0:
	test_loaded_strongest_on_strongest_component = "True"
else:
	test_loaded_strongest_on_strongest_component = "False"

t = open(output_file + 'analytics.txt','w')
t.write(str(test_tau) + '\t' + test_loaded_strongest_on_strongest_component + '\t' + str(corrz[0]) + '\n')

