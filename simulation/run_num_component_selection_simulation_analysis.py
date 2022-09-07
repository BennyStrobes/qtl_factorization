import numpy as np 
import os
import sys
import pdb
import surge_inference
import pandas as pd
from sklearn.linear_model import LinearRegression


# Generate genotype vector for a particular snp
def get_genotype_vec(num_samples, af):
	allele1 = np.random.binomial(n=1, p=af,size=num_samples)
	allele2 = np.random.binomial(n=1, p=af,size=num_samples)
	genotype_vec = allele1 + allele2
	return genotype_vec

# Generate gene expression data (Y), and genotype data (G)
def generate_eqtl_factorization_data_with_no_sample_repeat(num_samples, num_tests, num_components, t_statistic, missingness_fraction):
	# Simulate genotype
	G = np.zeros((num_samples, num_tests))
	G_raw = np.zeros((num_samples, num_tests))
	for test_num in range(num_tests):
		af = np.random.uniform(.1,.9)
		# Get genotype for this particular snp
		genotype_vec = get_genotype_vec(num_samples, af)
		G_raw[:, test_num] = genotype_vec
		# Standardized genotype
		G[:, test_num] = (genotype_vec - np.mean(genotype_vec))/np.std(genotype_vec)

	# Simulate U
	U = np.random.standard_normal(size=(num_samples, num_components))
	for component_num in range(num_components):
		U[:, component_num] = (U[:, component_num] - np.mean(U[:, component_num]))/np.std(U[:, component_num])
	
	# Simulate V
	V = np.random.normal(loc=0.0, scale=t_statistic, size=(num_components, num_tests))
	# Add missingness to V
	for test_num in range(num_tests):
		for component_num in range(num_components):
			V[component_num, test_num] = V[component_num, test_num]*np.random.binomial(n=1,p=missingness_fraction)

	# Get Expected expression value of each (sample, gene) pair
	predicted_mean = np.zeros((num_samples, num_tests))  + G*np.dot(U,V)
	# Simulate genotype fixed effect size
	betas = []
	Y = np.zeros((num_samples, num_tests))
	for test_num in range(num_tests):
		beta = np.random.normal(loc=0.0, scale=.1)
		predicted_mean[:, test_num] = predicted_mean[:, test_num] + G[:, test_num]*beta
		Y[:, test_num] = np.random.normal(predicted_mean[:,test_num])
		betas.append(beta)
	
	# Standardized Y
	for test_num in range(num_tests):
		Y[:, test_num] = (Y[:, test_num] - np.mean(Y[:, test_num]))/np.std(Y[:, test_num])
	return Y, G, np.asarray(betas), U, V



sim_data_dir = sys.argv[1]
eqtl_results_dir = sys.argv[2]
num_samples = int(sys.argv[3])
t_statistic = float(sys.argv[4])
missingness_fraction = float(sys.argv[5])


num_tests=1000
np.random.seed(2)

#num_samples = 250
#t_statistic = .5
#missingness_fraction = .3

simulated_ks = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]

num_reps = 10

rep_vec = []
num_sim_components_vec = []
num_learned_components_vec = []


for rep in range(num_reps):
	for simulated_k in simulated_ks:
		# Simulate data
		Y, G, betas, U_sim, V_sim = generate_eqtl_factorization_data_with_no_sample_repeat(num_samples, num_tests, simulated_k, t_statistic, missingness_fraction)

		# Run SURGE
		eqtl_vi = surge_inference.SURGE_VI(K=30, alpha=1e-3, beta=1e-3, ard_alpha=1e-3, ard_beta=1e-3, max_iter=6000, gamma_v=1.0, warmup_iterations=5, re_boolean=False, delta_elbo_threshold=.1, verbose=True, output_root=eqtl_results_dir)
		eqtl_vi.fit(G=G, G_fe=G, Y=Y, cov=np.ones((G.shape[0], 1)))

		# Extract relevent latent factors
		num_learned_comp = np.sum(eqtl_vi.factor_pve > 1e-5)

		print(num_learned_comp)
		rep_vec.append(rep)
		num_sim_components_vec.append(simulated_k)
		num_learned_components_vec.append(num_learned_comp)



df = pd.DataFrame({'repitition':rep_vec, 'number_simulated_components':num_sim_components_vec, 'number_learned_components':num_learned_components_vec})


df.to_csv(path_or_buf=eqtl_results_dir + 'number_of_learned_components_recovery_results_' + str(num_samples) + '_' + str(t_statistic) + '_' + str(missingness_fraction) + '.txt', sep='\t', index=False)

