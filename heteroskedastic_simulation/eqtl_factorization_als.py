import numpy as np 
import os
import sys
import pdb
import scipy.special as special
from sklearn import linear_model
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
import time
import sklearn.decomposition
import time



def outside_update_V_t(g_test, y_test, K, U, lasso_param):
	# Get U scaled by genotype for this test
	U_scaled = U*g_test[:,None]
		
	covariates = np.hstack((np.ones((U_scaled.shape[0],1)), np.transpose(np.asmatrix(g_test)), U_scaled))

	model = sm.OLS(y_test, covariates)
	alpha_param = np.zeros(covariates.shape[1]) + lasso_param
	alpha_param[0] = 0
	alpha_param[1] = 0
	fit = model.fit_regularized(method='elastic_net', alpha=alpha_param,L1_wt=0.0)

	return fit.params


def outside_update_U_n(g_sample, y_sample, K, V, lasso_param):
	# Get V scaled by genotype for this sample
	V_scaled = np.transpose(V)*g_sample[:,None]

	clf = linear_model.ElasticNet(alpha=lasso_param, l1_ratio=1.0, positive=True, fit_intercept=False)
	clf.fit(V_scaled, y_sample)

	return clf.coef_


class EQTL_FACTORIZATION_ALS(object):
	def __init__(self, K=5, max_iter=300, lasso_param_u=.0001, lasso_param_v=.0001):
		self.max_iter = max_iter
		self.K = K
		self.lasso_param_u = lasso_param_u
		self.lasso_param_v = lasso_param_v
		self.iter = 0
	def fit(self, G, Y, z, cov):
		""" Fit the model.
			Args:
			G: A genotype matrix of floats with shape [num_samples, num_tests].
			Y: An expression matrix of floats with shape [num_samples, num_tests].
			z: groupings of length num_samples
		"""
		self.Y = Y
		self.G = G
		self.Z = z
		self.cov = cov
		self.initialize_variables()
		#self.update_elbo()
		# Loop through VI iterations
		for iter_num in range(self.max_iter):
			start_time = time.time()
			# Update parameter estimaters via ALS
			self.update_V()
			self.update_U()
			self.iter = self.iter + 1

			print('Variational Inference iteration: ' + str(iter_num))
			####################
			end_time = time.time()
			print(end_time-start_time)
			print('##############')
			print('##############')

	def update_V(self):
		###################
		# UPDATE V
		###################
		# Keep track of variables
		V_update_data = []

		for test_index in range(self.T):
			V_update_data.append(outside_update_V_t(self.G[:, test_index], self.Y[:, test_index], self.K, self.U, self.lasso_param_v))

		# Convert to array
		V_update_data = np.asarray(V_update_data).T

		self.intercept = V_update_data[0,:]
		self.V = V_update_data[1:,:]
	def update_U(self):
		Y_scaled = self.Y - self.intercept - (self.V[0,:]*self.G)

		U_update_data = []
		###################
		# UPDATE U
		###################
		# Don't parrallelize
		for sample_index in range(self.N):
			U_update_data.append(outside_update_U_n(self.G[sample_index, :], Y_scaled[sample_index, :], self.K, self.V[1:,:], self.lasso_param_u))

		# Convert to array
		U_update_data = np.asarray(U_update_data)
		#self.sample_intercept = U_update_data[:,0]
		self.U = U_update_data


	def initialize_variables(self):
		# Do standard variational inference
		self.N = self.Y.shape[0]
		self.T = self.Y.shape[1]
		self.num_cov = self.cov.shape[1]
		self.U = np.random.random(size=(self.N, self.K))
		self.V = np.zeros((self.K, self.T))
		self.C = np.zeros((self.num_cov, self.T))
