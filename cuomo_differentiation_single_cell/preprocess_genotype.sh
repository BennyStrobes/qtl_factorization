#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=lrgmem
#SBATCH --mem=100GB

#SBATCH --nodes=1



genotype_dir="$1"
processed_genotype_dir="$2"


python preprocess_genotype.py $genotype_dir $processed_genotype_dir