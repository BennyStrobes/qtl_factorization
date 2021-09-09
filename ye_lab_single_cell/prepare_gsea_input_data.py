import numpy as np 
import os
import sys
import pdb







model_loadings_file = sys.argv[1]
expr_file = sys.argv[2]
gene_names_file = sys.argv[3]
component_num = int(sys.argv[4])
component_output_stem = sys.argv[5]


loadings = np.loadtxt(model_loadings_file)
expr = np.loadtxt(expr_file)
gene_names = np.loadtxt(gene_names_file, dtype=str, delimiter='\t')
gene_symbols = gene_names[:,0]

loading_vec = loadings[:, component_num]
corrz = []

num_genes = gene_names.shape[0]


for gene_num in range(num_genes):
	corrz.append(np.abs(np.corrcoef(expr[:, gene_num], loading_vec)[0,1]))

corrz = np.asarray(corrz)

ordered_gene_indices = np.argsort((-corrz))

num_test_genes = 50
test_genes = []

pdb.set_trace()

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



