import numpy as np 
import os
import sys
import pdb
import eqtl_factorization_ard_no_alpha
from sklearn.linear_model import LinearRegression


def get_genotype_vec(num_samples):
	maf = np.random.uniform(.05,.5)
	allele1 = np.random.binomial(n=1, p=maf,size=num_samples)
	allele2 = np.random.binomial(n=1, p=maf,size=num_samples)
	genotype_vec = allele1 + allele2
	return genotype_vec

def generate_null_eqtl_factorization_data_with_no_sample_repeat_by_variable_variance_ratios(num_samples, num_tests):
	# Simulate genotype
	G = np.zeros((num_samples, num_tests))
	G_raw = np.zeros((num_samples, num_tests))
	mafs = []
	for test_num in range(num_tests):
		# Get genotype for this particular snp
		genotype_vec = get_genotype_vec(num_samples)
		if np.mean(genotype_vec) == 1.0:
			genotype_vec = get_genotype_vec(num_samples)
			if np.mean(genotype_vec) == 1.0:
				genotype_vec = get_genotype_vec(num_samples)
				if np.mean(genotype_vec) == 1.0:
					genotype_vec = get_genotype_vec(num_samples)
		af = np.sum(genotype_vec)/(2.0*len(genotype_vec))
		if af > .5:
			maf = 1.0 - af
		else:
			maf = af
		mafs.append(maf)
		G_raw[:, test_num] = genotype_vec
		G[:, test_num] = (genotype_vec - np.mean(genotype_vec))/np.std(genotype_vec)
	mafs = np.asarray(mafs)
	covs = np.random.randn(num_samples,5)
	C = np.random.randn(5,num_tests)
	predicted_mean = np.zeros((num_samples, num_tests))
	used = False
	betas = []
	versions = []
	Y = np.zeros((num_samples, num_tests))
	for test_num in range(num_tests):
		if np.random.uniform() < .04:
			beta = np.random.normal(loc=0.0, scale=.1)
			beta = 0
			predicted_mean[:, test_num] = predicted_mean[:, test_num] + G[:, test_num]*beta
			versions.append(0)
			indices0 = G_raw[:,test_num] == 0.0
			indices1 = G_raw[:,test_num] == 1.0
			indices2 = G_raw[:,test_num] == 2.0
			Y[indices0,test_num] = np.random.normal(predicted_mean[indices0,test_num], scale=.1)
			Y[indices1,test_num] = np.random.normal(predicted_mean[indices1,test_num], scale=1)
			Y[indices2,test_num] = np.random.normal(predicted_mean[indices2,test_num], scale=5)
		elif np.random.uniform() < .04:
			beta = np.random.normal(loc=0.0, scale=.1)
			beta = 0
			predicted_mean[:, test_num] = predicted_mean[:, test_num] + G[:, test_num]*beta
			versions.append(1)
			indices0 = G_raw[:,test_num] == 0.0
			indices1 = G_raw[:,test_num] == 1.0
			indices2 = G_raw[:,test_num] == 2.0
			Y[indices0,test_num] = np.random.normal(predicted_mean[indices0,test_num], scale=5)
			Y[indices1,test_num] = np.random.normal(predicted_mean[indices1,test_num], scale=1)
			Y[indices2,test_num] = np.random.normal(predicted_mean[indices2,test_num], scale=.1)
		else:
			beta = np.random.normal(loc=0.0, scale=.1)
			beta = 0
			predicted_mean[:, test_num] = predicted_mean[:, test_num] + G[:, test_num]*beta
			Y[:, test_num] = np.random.normal(predicted_mean[:,test_num])
			versions.append(2)
		betas.append(beta)
	for test_num in range(num_tests):
		Y[:, test_num] = (Y[:, test_num] - np.mean(Y[:, test_num]))/np.std(Y[:, test_num])
	return Y, G, np.asarray(betas), np.asarray(versions), mafs

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
	return np.asarray(resids), np.asarray(learned_betas), np.asarray(corrz)

####################
# Command line args
####################
num_samples = int(sys.argv[1])
num_tests = int(sys.argv[2])



Y, G, betas, versions, mafs = generate_null_eqtl_factorization_data_with_no_sample_repeat_by_variable_variance_ratios(num_samples, num_tests)
Y_resid, beta_estimated, corrz = get_residuals(Y, G)
print(sum(versions==0))
print(sum(versions==1))

print(np.sort(np.abs(corrz)))
pdb.set_trace()
# Dummy variable
Z = np.ones(num_samples)
cov = np.ones((num_samples, 1))


eqtl_vi = eqtl_factorization_ard_no_alpha.EQTL_FACTORIZATION_VI(K=10, alpha=1e-3, beta=1e-3, a=1, b=1, max_iter=200, gamma_v=1)
eqtl_vi.fit(G=G, Y=Y, z=Z, cov=cov)
pdb.set_trace()