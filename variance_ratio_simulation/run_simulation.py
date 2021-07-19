import numpy as np 
import os
import sys
import pdb
import eqtl_factorization_ard_no_alpha

def get_genotype_vec(num_samples):
	maf = np.random.uniform(.1,.9)
	allele1 = np.random.binomial(n=1, p=maf,size=num_samples)
	allele2 = np.random.binomial(n=1, p=maf,size=num_samples)
	genotype_vec = allele1 + allele2
	return genotype_vec

def generate_null_eqtl_factorization_data_with_no_sample_repeat_by_variable_variance_ratios(num_samples, num_tests):
	# Simulate genotype
	G = np.zeros((num_samples, num_tests))
	for test_num in range(num_tests):
		# Get genotype for this particular snp
		genotype_vec = get_genotype_vec(num_samples)
		if np.mean(genotype_vec) == 1.0:
			genotype_vec = get_genotype_vec(num_samples)
			if np.mean(genotype_vec) == 1.0:
				genotype_vec = get_genotype_vec(num_samples)
				if np.mean(genotype_vec) == 1.0:
					genotype_vec = get_genotype_vec(num_samples)
		G[:, test_num] = (genotype_vec - np.mean(genotype_vec))/np.std(genotype_vec)
	
	covs = np.random.randn(num_samples,5)
	C = np.random.randn(5,num_tests)
	predicted_mean = np.dot(covs,C)
	for test_num in range(num_tests):
		beta = np.random.randn()
		predicted_mean[:, test_num] = predicted_mean[:, test_num] + G[:, test_num]*beta
	Y = np.random.normal(predicted_mean)
	for test_num in range(num_tests):
		Y[:, test_num] = (Y[:, test_num] - np.mean(Y[:, test_num]))/np.std(Y[:, test_num])
	return Y, G

####################
# Command line args
####################
num_samples = int(sys.argv[1])
num_tests = int(sys.argv[2])



Y, G = generate_null_eqtl_factorization_data_with_no_sample_repeat_by_variable_variance_ratios(num_samples, num_tests)
G[:,-1] = G[:,-1]*1000
# Dummy variable
Z = np.ones(num_samples)
cov = np.ones((num_samples, 1))


eqtl_vi = eqtl_factorization_ard_no_alpha.EQTL_FACTORIZATION_VI(K=10, alpha=1e-3, beta=1e-3, a=1, b=1, max_iter=100, gamma_v=1)
eqtl_vi.fit(G=G, Y=Y, z=Z, cov=cov)
pdb.set_trace()