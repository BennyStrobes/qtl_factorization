import numpy as np 
import os
import sys
import pdb





#####################
# commannd line args
#####################
qtl_covariate_file = sys.argv[1]
qtl_lf_file = sys.argv[2]
qtl_new_covariate_file = sys.argv[3]


cov1 = np.loadtxt(qtl_covariate_file)
cov2 = np.loadtxt(qtl_lf_file)

merged_cov = np.hstack((cov1,cov2))


np.savetxt(qtl_new_covariate_file, merged_cov, fmt="%s", delimiter='\t')