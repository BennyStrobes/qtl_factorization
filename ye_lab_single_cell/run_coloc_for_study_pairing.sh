#!/bin/bash -l

#SBATCH
#SBATCH --time=20:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1


sig_eqtl_file="$1"
all_eqtl_file="$2"
gwas_study_file_root="$3"
output_root="$4"
eqtl_study_name="$5"
gwas_study_name="$6"
eqtl_data_version="$7"

echo "preparing coloc data"
if false; then
python prepare_coloc_data_for_study_pairing.py $sig_eqtl_file $all_eqtl_file $gwas_study_file_root $output_root $eqtl_data_version
fi

echo "running coloc"
module load R/3.5.1
Rscript run_coloc_for_study_pairing.R $output_root $eqtl_study_name $gwas_study_name
