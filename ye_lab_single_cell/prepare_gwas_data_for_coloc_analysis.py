import numpy as np 
import os
import sys
import pdb
import gzip

def reprint_bentham_2015_sle_study(input_study_file, output_study_file_root):
	for chrom_num in range(1,23):
		chrom_string = str(chrom_num)
		f = gzip.open(input_study_file)
		t = open(output_study_file_root + chrom_string + '.txt','w')
		t.write('variant_id\tbeta\tvar_beta\tp_value\tN\n')
		head_count = 0
		for line in f:
			line = line.decode('UTF-8').rstrip()
			data = line.split('\t')
			if len(data) != 11:
				print('assumption eroror')
				pdb.set_trace()
			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue
			line_chrom_num = data[0]
			if line_chrom_num != chrom_string:
				continue
			# Extract relevant fields
			variant_id = data[0] + ':' + data[1] + ':' + data[3] + ':' + data[4]
			beta = data[6]
			se_beta = float(data[7])
			if se_beta == 0.0:
				continue
			var_beta = str(np.square(se_beta))
			pvalue = data[5]
			N_eff = '23210'
			# Print to output file
			t.write(variant_id + '\t' + beta + '\t' + var_beta + '\t' + pvalue + '\t' + N_eff + '\n')
		f.close()
		t.close()

def reprint_ukbb_study(input_study_file, output_study_file_root):
	for chrom_num in range(1,23):
		chrom_string = str(chrom_num)
		f = gzip.open(input_study_file)
		t = open(output_study_file_root + chrom_string + '.txt','w')
		t.write('variant_id\tbeta\tvar_beta\tp_value\tN\n')
		head_count = 0
		for line in f:
			line = line.decode('UTF-8').rstrip()
			data = line.split('\t')
			# Skip header
			if head_count == 0:
				head_count = head_count + 1
				continue
			# Simple error checking
			if len(data) != 12:
				print('assumption erorroror')
				pdb.set_trace()
			# Skip lines not on the current chromosome
			line_chrom_num = data[1]
			if line_chrom_num != chrom_string:
				continue
			# Extract relevant fields
			variant_id = data[1] + ':' + data[2] + ':' + data[3] + ':' + data[4]
			beta = data[7]
			se_beta = float(data[8])
			var_beta = str(np.square(se_beta))
			pvalue = data[9]
			N_eff = data[10]
			# Print to output file
			t.write(variant_id + '\t' + beta + '\t' + var_beta + '\t' + pvalue + '\t' + N_eff + '\n')
		f.close()
		t.close()


coloc_input_dir = sys.argv[1]
processed_gwas_studies_file = sys.argv[2]
coloc_dir = sys.argv[3]


# Open processed gwas studies file and print header
t_meta = open(processed_gwas_studies_file,'w')
#t_meta.write('study_name\tsumstats_path_root\n')


################
# UKBB Studies
################
# blood monocyte count
input_study_file = coloc_input_dir + 'blood_MONOCYTE_COUNT.sumstats.gz'
study_name = 'ukbb_blood_monocyte_count'
output_study_file_root = coloc_dir + study_name + '_coloc_sumstats_'
t_meta.write(study_name + '\t' + output_study_file_root + '\n')
#reprint_ukbb_study(input_study_file, output_study_file_root)


# blood lymphocyte count
input_study_file = coloc_input_dir + 'blood_LYMPHOCYTE_COUNT.sumstats.gz'
study_name = 'ukbb_blood_lymphocyte_count'
output_study_file_root = coloc_dir + study_name + '_coloc_sumstats_'
t_meta.write(study_name + '\t' + output_study_file_root + '\n')
#reprint_ukbb_study(input_study_file, output_study_file_root)

# eczema
input_study_file = coloc_input_dir + 'disease_ALLERGY_ECZEMA_DIAGNOSED.sumstats.gz'
study_name = 'ukbb_eczema'
output_study_file_root = coloc_dir + study_name + '_coloc_sumstats_'
t_meta.write(study_name + '\t' + output_study_file_root + '\n')
#reprint_ukbb_study(input_study_file, output_study_file_root)

# bmi
input_study_file = coloc_input_dir + 'body_BMIz.sumstats.gz'
study_name = 'ukbb_bmi'
output_study_file_root = coloc_dir + study_name + '_coloc_sumstats_'
t_meta.write(study_name + '\t' + output_study_file_root + '\n')
#reprint_ukbb_study(input_study_file, output_study_file_root)

input_study_file = coloc_input_dir + 'blood_EOSINOPHIL_COUNT.sumstats.gz'
study_name = 'ukbb_blood_eosinophil_count'
output_study_file_root = coloc_dir + study_name + '_coloc_sumstats_'
t_meta.write(study_name + '\t' + output_study_file_root + '\n')
#reprint_ukbb_study(input_study_file, output_study_file_root)

input_study_file = coloc_input_dir + 'blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT.sumstats.gz'
study_name = 'ukbb_blood_high_light_scatter_reticulotye_count'
output_study_file_root = coloc_dir + study_name + '_coloc_sumstats_'
t_meta.write(study_name + '\t' + output_study_file_root + '\n')
#reprint_ukbb_study(input_study_file, output_study_file_root)

input_study_file = coloc_input_dir + 'blood_MEAN_CORPUSCULAR_HEMOGLOBIN.sumstats.gz'
study_name = 'ukbb_blood_mean_corpuscular_hemoglobin'
output_study_file_root = coloc_dir + study_name + '_coloc_sumstats_'
t_meta.write(study_name + '\t' + output_study_file_root + '\n')
#reprint_ukbb_study(input_study_file, output_study_file_root)

input_study_file = coloc_input_dir + 'blood_MEAN_PLATELET_VOL.sumstats.gz'
study_name = 'ukbb_blood_platelet_vol'
output_study_file_root = coloc_dir + study_name + '_coloc_sumstats_'
t_meta.write(study_name + '\t' + output_study_file_root + '\n')
#reprint_ukbb_study(input_study_file, output_study_file_root)

input_study_file = coloc_input_dir + 'blood_PLATELET_COUNT.sumstats.gz'
study_name = 'ukbb_blood_platelet_count'
output_study_file_root = coloc_dir + study_name + '_coloc_sumstats_'
t_meta.write(study_name + '\t' + output_study_file_root + '\n')
#reprint_ukbb_study(input_study_file, output_study_file_root)

input_study_file = coloc_input_dir + 'blood_RED_COUNT.sumstats.gz'
study_name = 'ukbb_blood_red_count'
output_study_file_root = coloc_dir + study_name + '_coloc_sumstats_'
t_meta.write(study_name + '\t' + output_study_file_root + '\n')
#reprint_ukbb_study(input_study_file, output_study_file_root)

input_study_file = coloc_input_dir + 'blood_WHITE_COUNT.sumstats.gz'
study_name = 'ukbb_blood_white_count'
output_study_file_root = coloc_dir + study_name + '_coloc_sumstats_'
t_meta.write(study_name + '\t' + output_study_file_root + '\n')
#reprint_ukbb_study(input_study_file, output_study_file_root)

input_study_file = coloc_input_dir + 'bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv.gz'
study_name = 'sle'
output_study_file_root = coloc_dir + study_name + '_coloc_sumstats_'
t_meta.write(study_name + '\t' + output_study_file_root + '\n')
reprint_bentham_2015_sle_study(input_study_file, output_study_file_root)


t_meta.close()
