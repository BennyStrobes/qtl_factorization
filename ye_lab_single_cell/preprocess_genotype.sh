#!/bin/bash -l

#SBATCH
#SBATCH --time=40:00:00
#SBATCH --mem=4GB
#SBATCH --partition=lrgmem
#SBATCH --nodes=1


genotype_data_dir="$1"
processed_genotype_dir="$2"


bcftools merge $genotype_data_dir"all_clues.processed.vcf.gz" $genotype_data_dir"immvar.processed.vcf.gz" -Oz -o $processed_genotype_dir"clues_immvar_merged.vcf.gz"