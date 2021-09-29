import numpy as np 
import os
import sys
import pdb




surge_latent_factors_file = sys.argv[1]
factor_pve_file = sys.argv[2]
surge_latent_factor_output_file = sys.argv[3]


U = np.loadtxt(surge_latent_factors_file, dtype=str, delimiter='\t')

factor_pves = np.loadtxt(factor_pve_file, delimiter='\t')

ordered_factors = np.argsort(-factor_pves)


U_reordered = U[:, ordered_factors]

np.savetxt(surge_latent_factor_output_file, U_reordered, fmt="%s", delimiter='\t')