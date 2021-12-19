import numpy as np 
import os
import sys
import pdb



def extract_enrichment_and_tau_from_bottom_line_of_sldsc_results_file(per_cell_trait_file):
	aa = np.loadtxt(per_cell_trait_file,dtype=str, delimiter='\t')
	pvalue = float(aa[-1,6])
	neg_log10_pvalue = str(-np.log10(pvalue))
	enrichment = aa[-1,4]
	tau = aa[-1,-3]
	enrichment_std_err = aa[-1,5]
	tau_std_err = aa[-1,-2]
	return enrichment, tau, neg_log10_pvalue, enrichment_std_err, tau_std_err



def meta_analysis(effects, se, method='random', weights=None):
	# From Omer Weissbrod
	assert method in ['fixed', 'random']
	d = effects
	variances = se**2

	#compute random-effects variance tau2
	vwts = 1.0 / variances
	fixedsumm = vwts.dot(d) / vwts.sum()
	Q = np.sum(((d - fixedsumm)**2) / variances)
	df = len(d)-1
	tau2 = np.maximum(0, (Q-df) / (vwts.sum() - vwts.dot(vwts) / vwts.sum()))

	#defing weights
	if weights is None:
		if method == 'fixed':
			wt = 1.0 / variances
		else:
			wt = 1.0 / (variances + tau2)
	else:
		wt = weights

	#compute summtest
	summ = wt.dot(d) / wt.sum()
	if method == 'fixed':
		varsum = np.sum(wt*wt*variances) / (np.sum(wt)**2)
	else:
		varsum = np.sum(wt*wt*(variances+tau2)) / (np.sum(wt)**2)
	###summtest = summ / np.sqrt(varsum)

	summary=summ
	se_summary=np.sqrt(varsum)

	return summary, se_summary


def sldsc_meta_analysis(meta_analysis_name, trait_names, t, sample_names_file, per_cell_sldsc_results_dir):
	f = open(sample_names_file)
	for line in f:
		sample_name = line.rstrip()
		sample_measured_in_all_traits = True
		component_num = str(int(sample_name.split(':')[0].split('ent_')[1]) + 1)
		component_position = sample_name.split(':')[1]
		enrichments = []
		taus = []
		enrichment_ses = []
		tau_ses = []
		for trait_name in trait_names:
			per_cell_trait_file = per_cell_sldsc_results_dir + 'sample_specific_eqtl_effect_sizes_' + trait_name + '_' + sample_name + '.results'
			if os.path.exists(per_cell_trait_file) == False:
				sample_measured_in_all_traits = False
			else:
				enrichment, tau, neg_log10_pvalue, enrichment_std_err, tau_std_err = extract_enrichment_and_tau_from_bottom_line_of_sldsc_results_file(per_cell_trait_file)
				enrichments.append(float(enrichment))
				taus.append(float(tau))
				enrichment_ses.append(float(enrichment_std_err))
				tau_ses.append(float(tau_std_err))
		enrichments = np.asarray(enrichments)
		taus = np.asarray(taus)
		enrichment_ses = np.asarray(enrichment_ses)
		tau_ses = np.asarray(tau_ses)
		# Run meta analysis
		if sample_measured_in_all_traits == False:
			t.write(sample_name + '\t' + component_num + '\t' + component_position + '\t' + meta_analysis_name + '\tNA\tNA\tNA\tNA\tNA\n')
		else:
			ma_enrichment, ma_enrichment_std_error = meta_analysis(enrichments, enrichment_ses, method='random')
			ma_tau, ma_tau_std_error = meta_analysis(taus, tau_ses, method='random')
			t.write(sample_name + '\t' + component_num + '\t' + component_position + '\t' + meta_analysis_name + '\t' + str(ma_enrichment) + '\t' + str(ma_enrichment_std_error) + '\t' + str(ma_tau) + '\t' + str(ma_tau_std_error) + '\t' + 'NA' + '\n')
	f.close()
	return t




sample_names_file = sys.argv[1]
per_cell_sldsc_results_dir = sys.argv[2]

# Number of snps
M = 5961159

trait_names = ['ukbb_blood_monocyte_count', 'ukbb_blood_lymphocyte_count', 'ukbb_bmi', 'ukbb_eczema', 'ukbb_blood_eosinophil_count', 'ukbb_blood_high_light_scatter_reticulotye_count', 'ukbb_blood_mean_corpuscular_hemoglobin', 'ukbb_blood_platelet_vol', 'ukbb_blood_platelet_count', 'ukbb_blood_red_count', 'ukbb_blood_white_count', 'ukbb_height', 'ukbb_T2D', 'Celiac', 'Crohns', 'Ulcerative_Colitis', 'Rheumatoid_Arthritis', 'Lupus', 'IBD', 'Multiple_sclerosis', 'PBC', 'CAD', 'Bipolar', 'Alzheimer', 'Schizophrenia']

print('start')
output_file = per_cell_sldsc_results_dir + 'component_gridspace_sldsc_results.txt'

t = open(output_file,'w')
t.write('sample_name\tcomponent_num\tcomponent_position\ttrait_name\tenrichment\tenrichment_std_err\ttau_star\ttau_star_std_err\tneg_log10_pvalue\n')

f = open(sample_names_file)
count = 0
for line in f:
	count = count + 1
	sample_name = line.rstrip()
	component_num = str(int(sample_name.split(':')[0].split('ent_')[1]) + 1)
	component_position = sample_name.split(':')[1]
	for trait_name in trait_names:
		per_cell_trait_file = per_cell_sldsc_results_dir + 'sample_specific_eqtl_effect_sizes_' + trait_name + '_' + sample_name + '.results'
		if os.path.exists(per_cell_trait_file) == False:
			t.write(sample_name + '\t' + component_num + '\t' + component_position + '\t' + trait_name + '\tNA\tNA\tNA\tNA\tNA\n')
		else:
			enrichment, tau, neg_log10_pvalue, enrichment_std_err, tau_std_err = extract_enrichment_and_tau_from_bottom_line_of_sldsc_results_file(per_cell_trait_file)
			t.write(sample_name + '\t' + component_num + '\t' + component_position + '\t' + trait_name + '\t' + enrichment + '\t' + enrichment_std_err + '\t' + tau + '\t' + tau_std_err + '\t' + neg_log10_pvalue + '\n')

f.close()









meta_analysis_name = 'Blood_meta'
trait_names = ['ukbb_blood_eosinophil_count', 'ukbb_blood_platelet_count', 'ukbb_blood_platelet_vol', 'ukbb_blood_red_count', 'ukbb_blood_white_count']
t = sldsc_meta_analysis(meta_analysis_name, trait_names, t, sample_names_file, per_cell_sldsc_results_dir)


meta_analysis_name = 'Immune_meta'
trait_names = ['Celiac', 'Crohns', 'ukbb_eczema', 'Lupus', 'Rheumatoid_Arthritis', 'Ulcerative_Colitis']
meta_analysis_file = per_cell_sldsc_results_dir + 'per_cell_sldsc_' + meta_analysis_name + '_meta_analysis_results.txt'
t = sldsc_meta_analysis(meta_analysis_name, trait_names, t, sample_names_file, per_cell_sldsc_results_dir)


meta_analysis_name = 'Non_blood_immune_meta'
trait_names = ['ukbb_bmi', 'CAD', 'ukbb_height', 'Schizophrenia', 'ukbb_T2D']
meta_analysis_file = per_cell_sldsc_results_dir + 'per_cell_sldsc_' + meta_analysis_name + '_meta_analysis_results.txt'
t = sldsc_meta_analysis(meta_analysis_name, trait_names, t, sample_names_file, per_cell_sldsc_results_dir)
t.close()
