#!/bin/bash -l

#SBATCH
#SBATCH --time=40:00:00
#SBATCH --mem=4GB
#SBATCH --partition=lrgmem
#SBATCH --nodes=1


genotype_data_dir="$1"
processed_genotype_dir="$2"
cell_covariates_file="$3"

#########################
# Merge genotype vcfs from two different studies into 1
########################
merged_vcf=$processed_genotype_dir"clues_immvar_merged.vcf.gz"
if false; then
bcftools merge $genotype_data_dir"all_clues.processed.vcf.gz" $genotype_data_dir"immvar.processed.vcf.gz" -Oz -o $merged_vcf
fi

#########################
# Extract list of individuals we have sc-rna-seq for
########################
sc_rna_seq_individual_file=$processed_genotype_dir"sc_rna_seq_individual_list.txt"
if false; then
python generate_list_of_individuals_with_rna_seq.py $cell_covariates_file $sc_rna_seq_individual_file
fi

#########################
# Filter vcf file to individuals we have sc-rna-seq for
########################
donor_filtered_merged_vcf=$processed_genotype_dir"clues_immvar_donor_filter_merged.vcf.gz"
vcftools --gzvcf $merged_vcf --keep $sc_rna_seq_individual_file --recode --stdout | gzip -c > $donor_filtered_merged_vcf


#########################
# Filter vcf file to sites with maf and no missing
########################
donor_site_filtered_merged_vcf=$processed_genotype_dir"clues_immvar_donor_site_filter_merged.vcf.gz"
vcftools --gzvcf $donor_filtered_merged_vcf --remove-filtered-all --max-missing 1 --maf .1 --recode --stdout | gzip -c > $donor_site_filtered_merged_vcf


