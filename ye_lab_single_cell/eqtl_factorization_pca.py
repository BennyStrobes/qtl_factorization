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
		# Run model
		divy = self.Y/self.G
		_pca = PCA(n_components=self.K, svd_solver='arpack')
		svd_loadings = _pca.fit_transform(divy)
		np.savetxt(self.output_root + 'temper_U_S.txt', svd_loadings, fmt="%s", delimiter='\t')
