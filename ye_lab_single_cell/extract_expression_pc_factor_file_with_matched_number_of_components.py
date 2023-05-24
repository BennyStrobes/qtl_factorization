import numpy as np 
import os
import sys
import pdb







surge_latent_factor_file = sys.argv[1]  # used to get number of components
qtl_covariate_file = sys.argv[2]  # Important: assuming first couple of columns of qtl_covariate_file are expression PCs
qtl_expression_pc_file = sys.argv[3]
expression_pc_factor_file = sys.argv[4] # output file
new_covariate_file = sys.argv[5]

# Get number of components
tmp = np.loadtxt(surge_latent_factor_file)
n_comp = tmp.shape[1]

# Load in expression_pcs
expression_pcs = np.loadtxt(qtl_expression_pc_file)[:, :n_comp]

# Extract covariates
covariates = np.loadtxt(qtl_covariate_file)

# save to output
np.savetxt(expression_pc_factor_file, expression_pcs, fmt="%s", delimiter='\t')


# get new covariates
new_covariates = []
for col_iter in range(covariates.shape[1]):
	valid_col = True
	for col_iter2 in range(expression_pcs.shape[1]):
		if np.array_equal(covariates[:, col_iter], expression_pcs[:, col_iter2]):
			valid_col = False
	if valid_col:
		new_covariates.append(covariates[:, col_iter])

new_covariates = np.transpose(np.asarray(new_covariates))

# save to output
np.savetxt(new_covariate_file, new_covariates, fmt="%s", delimiter='\t')
