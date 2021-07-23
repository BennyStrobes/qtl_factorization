import numpy as np 
import os
import sys
import pdb


def extract_individuals_from_file(input_individual_file):
	f = open(input_individual_file)
	indis = []
	for line in f:
		line = line.rstrip()
		indis.append(line)
	return np.asarray(indis)

def extract_snp_names_from_file(input_snp_name_file, chrom_string):
	f = open(input_snp_name_file)
	snp_names = []
	for line in f:
		line = line.rstrip()
		data = line.split()
		if data[0] != chrom_string:
			print('assumption error!')
			pdb.set_trace()
		pos = int(data[1])
		if pos < 0:
			print('assumption error')
			pdb.set_trace()
		if pos > 250000000:
			print('assumption error')
			pdb.set_trace()
		snp_name = data[0] + ':' + data[1]
		snp_names.append(snp_name)
	return np.asarray(snp_names)

def generate_chromosome_specific_genotype(chrom_num, input_genotype_file, input_individual_file, input_snp_name_file, output_genotype_file):
	# Extract input data
	indis = extract_individuals_from_file(input_individual_file)
	snp_names = extract_snp_names_from_file(input_snp_name_file, str(chrom_num))
	genotype = np.transpose(np.loadtxt(input_genotype_file))
	genotype = genotype[1:,:]
	# Open output file 
	t = open(output_genotype_file, 'w')
	# Print header of output file
	t.write('snp_id\t' + '\t'.join(indis) + '\n')
	# Print each snp to outputfile
	num_snps = len(snp_names)
	for snp_num in range(num_snps):
		snp_name = snp_names[snp_num]
		genotype_vec = genotype[snp_num, :]
		t.write(snp_name + '\t' + '\t'.join(genotype_vec.astype(str)) + '\n')
	t.close()

def mean_inpute_genotype(output_genotype_file, output_mean_inputed_genotype_file):
	f = open(output_genotype_file)
	t = open(output_mean_inputed_genotype_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		variant_id = data[0]
		genotype = np.asarray(data[1:]).astype(float)
		observed_indices = np.where(genotype!=-1)[0]
		observed_mean = np.mean(genotype[observed_indices])
		genotype_inputed = np.copy(genotype)
		for i, val in enumerate(genotype):
			if val == -1:
				genotype_inputed[i] = observed_mean
		t.write(variant_id + '\t' + '\t'.join(genotype_inputed.astype(str)) +  '\n')
	t.close()
	f.close()

def remove_missing_from_genotype(output_genotype_file, output_missing_removed_genotype_file):
	f = open(output_genotype_file)
	t = open(output_missing_removed_genotype_file, 'w')
	head_count = 0
	for line in f:
		line = line.rstrip()
		data = line.split()
		if head_count == 0:
			head_count = head_count + 1
			t.write(line + '\n')
			continue
		variant_id = data[0]
		genotype = np.asarray(data[1:]).astype(float)
		missing_indices = np.where(genotype==-1)[0]
		if len(missing_indices) > 0:
			continue
		t.write(line + '\n')
	t.close()
	f.close()

def merge_genotype_across_chromosomes(genotype_input_stem, genotype_output):
	t = open(genotype_output, 'w')
	for chrom_num in range(1,23):
		chromosome_genotype_file = genotype_input_stem + str(chrom_num) + '.txt'
		f = open(chromosome_genotype_file)
		head_count = 0
		for line in f:
			line = line = line.rstrip()
			data = line.split()
			if head_count == 0:
				head_count = head_count + 1
				if chrom_num == 1:
					t.write(line + '\n')
					header = np.asarray(data)
				else:
					if np.array_equal(header, np.asarray(data)) == False:
						print('assumption error!')
						pdb.set_trace()
				continue
			t.write(line + '\n')
		f.close()
	t.close()

def compute_maf(genotype):
	af = np.sum(genotype)/(2.0*len(genotype))
	if af > .5:
		maf = 1.0 - af
	else:
		maf = af
	if maf > 0.5 or maf < 0.0:
		print('genotype assumption error')
		pdb.set_trace()
	return maf

#####################
# Command Line Args
#####################
genotype_dir = sys.argv[1]
pre_processed_data_dir = sys.argv[2]

# Generate chromosome-specific genotype files
for chrom_num in range(1,23):
	print(chrom_num)
	input_root = genotype_dir + 'chr' + str(chrom_num) + '/'
	input_genotype_file = input_root + 'combined.chr' + str(chrom_num) + '.common.012'
	input_individual_file = input_root + 'combined.chr' + str(chrom_num) + '.common.012.indv'
	input_snp_name_file = input_root + 'combined.chr' + str(chrom_num) + '.common.012.pos'
	output_genotype_file = pre_processed_data_dir + 'genotype_chr_' + str(chrom_num) + '.txt'
	generate_chromosome_specific_genotype(chrom_num, input_genotype_file, input_individual_file, input_snp_name_file, output_genotype_file)
	output_mean_inputed_genotype_file = pre_processed_data_dir + 'genotype_mean_inputed_chr_' + str(chrom_num) + '.txt'
	mean_inpute_genotype(output_genotype_file, output_mean_inputed_genotype_file)
	output_missing_removed_genotype_file = pre_processed_data_dir + 'genotype_missing_removed_chr_' + str(chrom_num) + '.txt'
	remove_missing_from_genotype(output_genotype_file, output_missing_removed_genotype_file)



# Generate genotype files across all chromosomes
genotype_input_stem = pre_processed_data_dir + 'genotype_chr_'
genotype_output = pre_processed_data_dir + 'genotype.txt'
merge_genotype_across_chromosomes(genotype_input_stem, genotype_output)


# Generate genotype files across all chromosomes
genotype_input_stem = pre_processed_data_dir + 'genotype_mean_inputed_chr_'
genotype_output = pre_processed_data_dir + 'genotype_mean_inputed.txt'
merge_genotype_across_chromosomes(genotype_input_stem, genotype_output)


# Generate genotype files across all chromosomes
genotype_input_stem = pre_processed_data_dir + 'genotype_missing_removed_chr_'
genotype_output = pre_processed_data_dir + 'genotype_missing_removed.txt'
merge_genotype_across_chromosomes(genotype_input_stem, genotype_output)