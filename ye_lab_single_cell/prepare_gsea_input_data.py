import numpy as np 
import os
import sys
import pdb
from sklearn.linear_model import LinearRegression
import scipy.stats






model_loadings_file = sys.argv[1]
model_pve_file = sys.argv[2]
expr_file = sys.argv[3]
gene_names_file = sys.argv[4]
component_num = int(sys.argv[5])
component_output_stem = sys.argv[6]


loadings = np.loadtxt(model_loadings_file)
expr = np.loadtxt(expr_file)
gene_names = np.loadtxt(gene_names_file, dtype=str, delimiter='\t')
gene_symbols = gene_names[:,0]

pve = np.loadtxt(model_pve_file)

ordered_components = np.argsort(-pve)
ordered_loadings = loadings[:, ordered_components]
loading_vec = ordered_loadings[:, component_num]
corrz = []
corrz2 = []
pvalz = []

num_genes = gene_names.shape[0]

temp_arr = []
for temp_component_num in range(ordered_loadings.shape[1]):
	if temp_component_num != component_num:
		temp_arr.append(temp_component_num)
temp_arr = np.asarray(temp_arr)

#cell_type_loadings = ordered_loadings[:,:3]
cell_type_loadings = ordered_loadings[:, temp_arr]


for gene_num in range(num_genes):
	y_vec = expr[:, gene_num]
	reg = LinearRegression(fit_intercept=True).fit(cell_type_loadings, np.transpose(np.asmatrix(y_vec)))
	pred_y = reg.predict(cell_type_loadings)[:,0]
	resid_y = y_vec - pred_y
	corrz.append(np.abs(np.corrcoef(resid_y, loading_vec)[0,1]))
	corrz2.append((np.corrcoef(resid_y, loading_vec)[0,1]))
	pvalz.append(scipy.stats.pearsonr(resid_y,loading_vec)[1])

corrz = np.asarray(corrz)
corrz2 = np.asarray(corrz2)
pvalz = np.asarray(pvalz)


ordered_gene_indices = np.argsort((-corrz))

print(component_num)
print(-np.sort(-corrz)[:5])

gene_correlations_file = component_output_stem + 'gene_correlations.txt'
t = open(gene_correlations_file, 'w')
t.write('gene_name\tcorrelation\tabsolute_correlation\tcorrelation_pvalue\n')
for gene_index in ordered_gene_indices:
	t.write(gene_symbols[gene_index] + '\t' + str(corrz2[gene_index]) + '\t' + str(corrz[gene_index]) + '\t' + str(pvalz[gene_index]) + '\n')
t.close()
print('done')


num_test_genes = 50
test_genes = []


for gene_index in ordered_gene_indices[:num_test_genes]:
	test_genes.append(gene_symbols[gene_index])

bgrd_genes = []
for gene_index in ordered_gene_indices[500:]:
	bgrd_genes.append(gene_symbols[gene_index])

test_genes = np.asarray(test_genes)
bgrd_genes = np.asarray(bgrd_genes)

test_names_file = component_output_stem + 'test_genes.txt'
bgrd_genes_file = component_output_stem + 'background_genes.txt'


np.savetxt(test_names_file, test_genes,fmt="%s", delimiter='\n')
np.savetxt(bgrd_genes_file, bgrd_genes,fmt="%s", delimiter='\n')



