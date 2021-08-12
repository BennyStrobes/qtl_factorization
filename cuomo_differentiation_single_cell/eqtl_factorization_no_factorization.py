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

def sigmoid_function(x):
	return 1.0/(1.0 + np.exp(-x))


class EQTL_FACTORIZATION(object):
	def __init__(self, output_root=''):
		# Output root (directory) to save intermediate results
		self.output_root = output_root
	def fit(self, G, G_fe, Y, z, cov):
		""" Fit the model.
			Args:
			G: A genotype matrix of floats with shape [num_samples, num_tests].
			G_fe: A genotype matrix (used for fixed effect of genotype) of floats with shape [num_samples, num_tests]. Will be the same as G unless in a particular permutation scheme
			Y: An expression matrix of floats with shape [num_samples, num_tests].
			z: groupings of length num_samples
			cov: A covariate matrix of floats with shape [num_samples, num_covariates]  ... we assume this contains an intercept term
		"""
		self.G = G
		self.G_fe = G_fe
		self.Y = Y
		self.z = np.asarray(z)
		self.cov = cov
		# Initialize variables
		print('Initialize variables')
		self.initialize_variables()

		# Loop through tests
		for test_num in range(self.T):
			print(test_num)
			y_vec = self.Y[:,test_num]
			g_vec = self.G_fe[:,test_num]
			# Save lmm input files to text
			np.savetxt('y.txt', y_vec, fmt="%s", delimiter='\t')
			np.savetxt('g.txt', g_vec, fmt="%s", delimiter='\t')
			np.savetxt('cov.txt', self.cov, fmt="%s", delimiter='\t')
			np.savetxt('z.txt', self.z, fmt="%s", delimiter='\t')
			# Run LMM in R
			os.system('Rscript lmm.R')
			# Re-import results
			temp_alpha_mu = np.loadtxt('alpha_mu.txt')
			temp_C_mu = np.loadtxt('C_mu.txt')
			temp_F_mu = np.loadtxt('F_mu.txt')
			re_var = np.loadtxt('re_variance.txt')
			residual_var = np.loadtxt('residual_variance.txt')
			# Save to data structure
			self.alpha_mu[:, test_num] = temp_alpha_mu
			self.C_mu[:, test_num] = temp_C_mu
			self.F_mu[test_num] = temp_F_mu
			self.tau[test_num] = 1.0/residual_var
			self.psi[test_num] = 1.0/re_var
	def initialize_variables(self):
		# Add model dimensions to object
		self.N = self.Y.shape[0]
		self.T = self.Y.shape[1]
		self.num_cov = self.cov.shape[1]

		# Random effects
		self.I = len(np.unique(self.z))

		# Random effects
		self.alpha_mu = np.zeros((self.I, self.T))

		# Random effects variances
		self.psi = np.ones(self.T)

		# Initialize C and F
		self.F_mu = np.zeros(self.T)
		self.C_mu = np.zeros((self.num_cov, self.T))

		# Variances
		self.tau = np.ones(self.T)
		self.print_diagnostic_data()
	def print_diagnostic_data(self):
		print(str(self.N) + ' samples detected')
		print(str(self.T) + ' tests detected')

