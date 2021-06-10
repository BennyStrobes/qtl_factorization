#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --partition=lrgmem
#SBATCH --nodes=1


genotype_data_dir="$1"
processed_genotype_dir="$2"
cell_covariates_file="$3"
visualize_processed_genotype_dir="$4"

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
if false; then
vcftools --gzvcf $merged_vcf --keep $sc_rna_seq_individual_file --recode --stdout | gzip -c > $donor_filtered_merged_vcf
fi

#########################
# Filter vcf file to sites with maf .1 and no missing
########################
donor_site_filtered_merged_vcf=$processed_genotype_dir"clues_immvar_donor_site_filter_merged.vcf.gz"
if false; then
vcftools --gzvcf $donor_filtered_merged_vcf --remove-filtered-all --max-missing 1 --maf .1 --recode --stdout | gzip -c > $donor_site_filtered_merged_vcf
fi

#########################
# Output dosage matrix for each chromosome
########################
if false; then
for chrom_num in $(seq 1 22); do 
	chromosome_dosage_prefix=$processed_genotype_dir"clues_immvar_chrom_"$chrom_num
	vcftools --gzvcf $donor_site_filtered_merged_vcf --chr $chrom_num --extract-FORMAT-info DS --out $chromosome_dosage_prefix
done
fi



#########################
# Compute genotype PCs with PLINK
# Referenced here: https://choishingwan.github.io/PRS-Tutorial/plink/
########################
echo "GENOTYPE PCS"
module load plink/1.90b6.4
donor_site_filtered_merged_plink_stem=$processed_genotype_dir"clues_immvar_donor_site_filter_merged_plink"

# Convert from to plink format
if false; then
plink --vcf $donor_site_filtered_merged_vcf --out $donor_site_filtered_merged_plink_stem
fi

# Filter to independent sites
donor_site_filtered_merged_independent_plink_stem=$processed_genotype_dir"clues_immvar_donor_site_filter_merged_independent_plink"
if false; then
plink --bfile $donor_site_filtered_merged_plink_stem --indep-pairwise 200 50 0.25 --out $donor_site_filtered_merged_independent_plink_stem
fi

# Compute Genotype PCs
if false; then
plink --bfile $donor_site_filtered_merged_plink_stem --extract $donor_site_filtered_merged_independent_plink_stem".prune.in" --pca 6 --out $donor_site_filtered_merged_independent_plink_stem
fi

# Tidy up genotype PCs (at indiviudal level and cell level)
genotype_pcs_file=$donor_site_filtered_merged_independent_plink_stem".eigenvec"
module load python/3.7.4-anaconda

python add_genotype_pcs_to_cell_covariates_file.py $cell_covariates_file $sc_rna_seq_individual_file $genotype_pcs_file $processed_genotype_dir



#########################
# Visualize processed genotypes
########################
if false; then
module load R/3.5.1
Rscript visualize_processed_genotypes.R $processed_genotype_dir $visualize_processed_genotype_dir
fi
