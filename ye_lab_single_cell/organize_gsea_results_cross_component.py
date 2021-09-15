import numpy as np 
import os
import pdb
import sys


def get_significant_pathways(gene_set_enrichment_output_stem, gene_set_name, components):
	pathway_pvalues = {}
	pathway_oddsratios = {}
	for component in components:
		file_name = gene_set_enrichment_output_stem + 'component_' + str(component) + '_' + gene_set_name
		f = open(file_name)
		count = 0
		for line in f:
			line = line.rstrip()
			if count < 4:
				count = count + 1
				continue
			data = line.split('\t')
			pathway_name = data[0]
			corrected_pvalue = float(data[7])
			if corrected_pvalue < .05:
				pathway_pvalues[pathway_name] = np.ones(len(components))
				pathway_oddsratios[pathway_name] = np.ones(len(components))
		f.close()
	for i, component in enumerate(components):
		file_name = gene_set_enrichment_output_stem + 'component_' + str(component) + '_' + gene_set_name
		f = open(file_name)
		count = 0
		for line in f:
			line = line.rstrip()
			if count < 4:
				count = count + 1
				continue
			data = line.split('\t')
			pathway_name = data[0]
			corrected_pvalue = float(data[7])
			oddsratio = data[5]
			if pathway_name not in pathway_pvalues:
				continue
			pathway_pvalues[pathway_name][i] = corrected_pvalue
			pathway_oddsratios[pathway_name][i] = oddsratio
	return pathway_pvalues, pathway_oddsratios




gene_set_enrichment_output_stem = sys.argv[1]
gene_set_name = sys.argv[2]

components = np.arange(3,10)


pathway_pvalues, pathway_oddsratios = get_significant_pathways(gene_set_enrichment_output_stem, gene_set_name, components)


pathways = sorted(pathway_pvalues.keys())



t = open(gene_set_enrichment_output_stem + 'cross_component_significant_enrichments_' + gene_set_name,'w')

t.write('pathway\tcomponent_num\tadjusted_pvalue\todds_ratio\n')

for pathway in pathways:
	for i,component in enumerate(components):
		t.write(pathway + '\t' + str(component) + '\t' + str(pathway_pvalues[pathway][i]) + '\t' + str(pathway_oddsratios[pathway][i]) + '\n')
t.close()


