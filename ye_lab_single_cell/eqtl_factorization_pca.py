import numpy as np 
import os
import sys
import pdb
import scipy.special as special
from sklearn.linear_model import LinearRegression
import time
import sklearn.decomposition
from joblib import Parallel, delayed
import multiprocessing
import time
import pandas as pd 
from pymer4.models import Lmer
from sklearn.decomposition import PCA
from ppca import PPCA

def run_pca(rat, K):
	scaled_rat = scale_allelic_ratios(rat)
	ppca = PPCA()
	ppca.fit(data=np.transpose(scaled_rat), d=K, verbose=True)
	U = ppca.C
	V = ppca.transform()
	return U

# scale each column
def scale_allelic_ratios(rat):
	scaled_rat = np.copy(rat)
	N, P = scaled_rat.shape
	# loop through columns
	for column_index in range(P):
		column_mean = np.nanmean(rat[:,column_index])
		column_std = np.nanstd(rat[:,column_index])
		for sample_index in range(N):
			if np.isnan(rat[sample_index, column_index]) == False:
				standardized_rat = (rat[sample_index, column_index] - column_mean)/column_std
				scaled_rat[sample_index, column_index] = standardized_rat
	return scaled_rat

def extract_residuals_from_standard_eqtl_model(Y, G, cov, z):
	num_tests = Y.shape[1]
	F_betas = []
	C_betas = []
	residuals = []
	for test_number in range(num_tests):
		y_vec = Y[:,test_number]
		g_vec = G[:,test_number]
		X = np.hstack((np.transpose(np.asmatrix(g_vec)), cov))
		#dd = {'y':y_vec, 'g':g_vec, 'z':z}
		#df = pd.DataFrame(dd)
		#model = Lmer('y ~ g  + (1|z)', data=df)
		#pdb.set_trace()
		reg = LinearRegression(fit_intercept=False).fit(X, np.transpose(np.asmatrix(y_vec)))
		F_betas.append(reg.coef_[0][0])
		C_betas.append(reg.coef_[0][1:])
		pred_y = reg.predict(X)
		test_resid = y_vec - pred_y[:,0]
		residuals.append(test_resid)
	residuals = np.transpose(np.asarray(residuals))
	return residuals

def extract_genotype_residuals(residuals, G):
	num_tests = G.shape[1]
	num_samples = G.shape[0]
	geno_residuals = []
	for test_number in range(num_tests):
		test_geno = G[:, test_number]
		rounded_genotypes = np.round(test_geno)
		unique_genotypes = np.unique(rounded_genotypes)
		if len(unique_genotypes) > 4:
			continue
		for unique_genotype in unique_genotypes:
			indices = np.where(rounded_genotypes == unique_genotype)[0]
			# Initialize geno residuals
			geno_residual = np.asarray([np.nan]*num_samples)

			geno_residual[indices] = residuals[indices, test_number]
			geno_residuals.append(geno_residual)
	geno_residuals = np.transpose(np.asarray(geno_residuals))
	return geno_residuals

class EQTL_FACTORIZATION_PCA(object):
	def __init__(self, K=25, output_root=''):
		self.K = K
		self.output_root = output_root
	def fit(self, G, Y, z, cov):
		""" Fit the model.
			Args:
			G: A genotype matrix of floats with shape [num_samples, num_tests].
			Y: An expression matrix of floats with shape [num_samples, num_tests].
			z: groupings of length num_samples
			cov: A covariate matrix of floats with shape [num_samples, num_covariates]  ... we assume this contains an intercept term
		"""
		self.G = G
		self.Y = Y
		self.z = np.asarray(z)
		self.cov = cov
		# Compute residuals from standard eqtl model
		residuals = extract_residuals_from_standard_eqtl_model(self.Y, self.G, self.cov, self.z)

		# Construct genotype residuals
		geno_residuals = extract_genotype_residuals(residuals, self.G)


		residual_loadings = run_pca(residuals, 10)

		geno_residual_loadings = run_pca(geno_residuals, 10)

		np.savetxt(self.output_root + 'temper_U_S_null.txt', residual_loadings, fmt="%s", delimiter='\t')
		np.savetxt(self.output_root + 'temper_U_S.txt', geno_residual_loadings, fmt="%s", delimiter='\t')
		# Run model
		#divy = self.Y/self.G
		#_pca = PCA(n_components=self.K, svd_solver='arpack')
		#svd_loadings = _pca.fit_transform(divy)
		#np.savetxt(self.output_root + 'temper_U_S.txt', svd_loadings, fmt="%s", delimiter='\t')
