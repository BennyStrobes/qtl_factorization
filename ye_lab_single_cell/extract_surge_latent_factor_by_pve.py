import numpy as np 
import os
import sys
import pdb





surge_results_stem = sys.argv[1]
surge_latent_factor_output_file = sys.argv[2]
perm_surge_latent_factor_factor_file = sys.argv[3]
surge_results_suffix = sys.argv[4]

# Real Surge loadings
surge_U_real_file_name = surge_results_stem + 'False' + surge_results_suffix + 'U_S.txt'
# Real Surge PVes
surge_pve_real_file_name = surge_results_stem + 'False' + surge_results_suffix +'factor_pve.txt'

# Perm Surge loadings
surge_U_perm_file_name = surge_results_stem + 'interaction_only' + surge_results_suffix +'U_S.txt'
# Perm Surge PVes
surge_pve_perm_file_name = surge_results_stem + 'interaction_only' + surge_results_suffix + 'factor_pve.txt'

# SURGE loadings
U_real = np.loadtxt(surge_U_real_file_name, dtype=str, delimiter='\t')
U_perm = np.loadtxt(surge_U_perm_file_name, dtype=str, delimiter='\t')

# SURGE PVEs
pve_real = np.loadtxt(surge_pve_real_file_name, delimiter='\t')
pve_perm = np.loadtxt(surge_pve_perm_file_name, delimiter='\t')


# Number of contexts to include
num_contexts = np.sum(pve_real >= 2e-5)
print(num_contexts)


# Factor re-ordering
ordered_real_factors = np.argsort(-pve_real)[:num_contexts]
ordered_perm_factors = np.argsort(-pve_perm)[:num_contexts]

# Re-order data 
pve_real_reorder = pve_real[ordered_real_factors]
pve_perm_reorder = pve_perm[ordered_perm_factors]

U_real_reorder = U_real[:, ordered_real_factors]
U_perm_reorder = U_perm[:, ordered_perm_factors]


# Save filtered loadings
np.savetxt(surge_latent_factor_output_file, U_real_reorder, fmt="%s", delimiter='\t')
np.savetxt(perm_surge_latent_factor_factor_file, U_perm_reorder, fmt="%s", delimiter='\t')

