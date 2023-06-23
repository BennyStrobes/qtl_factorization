import numpy as np 
import os
import sys
import pdb



def extract_surge_gene_trait_pairs(trait_names, num_components, pip_threshold, coloc_dir, eqtl_stem):
	hits = {}
	for component_iter in range(1, (num_components+1)):
		for trait_name in trait_names:
			filer = coloc_dir + eqtl_stem + '_' + str(component_iter) + '_interaction_' + trait_name + '_coloc_test_results.txt'
			f = open(filer)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				gene_name = data[0]
				pph4 = float(data[-1])
				if pph4 < pip_threshold:
					continue
				if gene_name + ':' + trait_name in hits:
					hits[gene_name + ':' + trait_name][0].append(pph4)
					hits[gene_name + ':' + trait_name][1].append(component_iter)
				else:
					hits[gene_name + ':' + trait_name] = ([pph4], [component_iter])
			f.close()
	return hits

def fill_in_dictionary(merged_trait_pairs, indexer, trait_names, num_components, coloc_dir, eqtl_stem):
	for component_iter in range(1, (num_components+1)):
		for trait_name in trait_names:
			filer = coloc_dir + eqtl_stem + '_' + str(component_iter) + '_interaction_' + trait_name + '_coloc_test_results.txt'
			f = open(filer)
			head_count = 0
			for line in f:
				line = line.rstrip()
				data = line.split('\t')
				if head_count == 0:
					head_count = head_count + 1
					continue
				gene_name = data[0]
				pph4 = float(data[-1])
				if gene_name + ':' + trait_name in merged_trait_pairs:
					merged_trait_pairs[gene_name + ':' + trait_name][indexer] = np.max([merged_trait_pairs[gene_name + ':' + trait_name][indexer], pph4])
			f.close()
	return merged_trait_pairs

gwas_studies_file = sys.argv[1]
coloc_dir = sys.argv[2]

pip_threshold = .95

trait_names = np.loadtxt(gwas_studies_file,dtype=str,delimiter='\t')[:-1,0]
num_components = 6

surge_gene_trait_pairs_discovery = extract_surge_gene_trait_pairs(trait_names, num_components, .95, coloc_dir, 'surge_latent_factor')
expression_pc_gene_trait_pairs_discovery = extract_surge_gene_trait_pairs(trait_names, num_components, .95, coloc_dir, 'expression_pc')

surge_gene_trait_pairs_validation = extract_surge_gene_trait_pairs(trait_names, num_components, .1, coloc_dir, 'surge_latent_factor')
expression_pc_gene_trait_pairs_validation = extract_surge_gene_trait_pairs(trait_names, num_components, .1, coloc_dir, 'expression_pc')


shared = {}
for ele in [*surge_gene_trait_pairs_discovery]:
	if ele in expression_pc_gene_trait_pairs_discovery:
		shared[ele] = 1

print('# discovered by SURGE: ' + str(len(surge_gene_trait_pairs_discovery)))
print('# discovered by expression pc: ' + str(len(expression_pc_gene_trait_pairs_discovery)))
print('# discovered by both: ' + str(len(shared)))


unique_to_surge = {}
for ele in [*surge_gene_trait_pairs_discovery]:
	if ele not in expression_pc_gene_trait_pairs_validation:
		unique_to_surge[ele] = surge_gene_trait_pairs_discovery[ele][1]

unique_to_expression_pc = {}
for ele in [*expression_pc_gene_trait_pairs_discovery]:
	if ele not in surge_gene_trait_pairs_validation:
		unique_to_expression_pc[ele] = 1


print('Unique to SURGE: ' + str(len(unique_to_surge)))
print('Unique to expression pc: ' + str(len(unique_to_expression_pc)))

counter = 0
unique_to_surge_non_pc_components = {}
for ele in [*unique_to_surge]:
	comps = unique_to_surge[ele]
	passer = False
	for comp in comps:
		if comp == 3 or comp == 5 or comp == 6:
			passer = True
	if passer:
		unique_to_surge_non_pc_components[ele] = 1
		if ele.split(':')[1].startswith('ukbb_blood'):
			counter = counter + 1
print('Unique to SURGE from non-pc correlated comps ' + str(len(unique_to_surge_non_pc_components)))


# Extract set of trait pairs found in either analysis
merged_trait_pairs = {}
for ele in [*surge_gene_trait_pairs_validation]:
	merged_trait_pairs[ele] = np.zeros(2)
for ele in [*expression_pc_gene_trait_pairs_validation]:
	merged_trait_pairs[ele] = np.zeros(2)

merged_trait_pairs = fill_in_dictionary(merged_trait_pairs, 0, trait_names, num_components, coloc_dir, 'surge_latent_factor')
merged_trait_pairs = fill_in_dictionary(merged_trait_pairs, 1, trait_names, num_components, coloc_dir, 'expression_pc')

output_file = coloc_dir + 'surge_vs_expression_pc_coloc_hits.txt'
t = open(output_file,'w')
t.write('gene_name\ttrait_name\tgene_trait_name\tmax_surge_pph4\tmax_expression_pc_pph4\tblood_trait\n')
for ele in [*merged_trait_pairs]:
	gene_name = ele.split(':')[0]
	trait_name = ele.split(':')[1]
	blood_trait = 'False'
	if trait_name.startswith('ukbb_blood'):
		blood_trait = 'True'
	t.write(gene_name + '\t' + trait_name + '\t' + ele + '\t' + str(merged_trait_pairs[ele][0]) + '\t' + str(merged_trait_pairs[ele][1]) + '\t' + blood_trait + '\n')
t.close()


