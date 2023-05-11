import numpy as np 
import os
import sys
import pdb







surge_latent_factor_file = sys.argv[1]  # used to get number of components
qtl_covariate_file = sys.argv[2]  # Important: assuming first couple of columns of qtl_covariate_file are expression PCs
expression_pc_factor_file = sys.argv[3] # output file
new_covariate_file = sys.argv[4]

# Get number of components
tmp = np.loadtxt(surge_latent_factor_file)
n_comp = tmp.shape[1]


# Extract covariates
covariates = np.loadtxt(qtl_covariate_file)

# Get top expression pcs
expression_pcs = covariates[:,:n_comp]
# save to output
np.savetxt(expression_pc_factor_file, expression_pcs, fmt="%s", delimiter='\t')


# get new covariates
new_covariates = covariates[:, n_comp:]
# save to output
np.savetxt(new_covariate_file, new_covariates, fmt="%s", delimiter='\t')
