import numpy as np 
import os
import sys
import pdb










surge_latent_factors_file = sys.argv[1]
factor_pve_file = sys.argv[2]
latent_factor_num = int(sys.argv[3])
surge_latent_factor_output_file = sys.argv[4]


U = np.loadtxt(surge_latent_factors_file, dtype=str, delimiter='\t')

factor_pves = np.loadtxt(factor_pve_file, delimiter='\t')

ordered_factors = np.argsort(-factor_pves)

temp_latent_factor_num = ordered_factors[(latent_factor_num-1)]  # Minus 1 to go from base 1 to base 0

U_slice = U[:, temp_latent_factor_num]

np.savetxt(surge_latent_factor_output_file, U_slice, fmt="%s", delimiter='\t')