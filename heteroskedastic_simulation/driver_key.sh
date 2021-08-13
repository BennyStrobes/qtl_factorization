#!/bin/bash -l

#SBATCH
#SBATCH --time=10:00:00
#SBATCH --mem=20GB
#SBATCH --partition=lrgmem
#SBATCH --nodes=1



output_root="/work-zfs/abattle4/bstrober/qtl_factorization/heteroskedastic_simulation/"
simulation_results_dir=$output_root"simulation_results/"
visualization_dir=$output_root"visualization/"


num_samples="2000"
num_tests="100"
version="heteroskedasticity"
af_array=(".1" ".2" ".3" ".4" ".5" ".6" ".7" ".8" ".9")
seed_array=("1" "2" "3" "4" "5")
for af in "${af_array[@]}"; do
	for seed in "${seed_array[@]}"; do
		output_file=$simulation_results_dir"results_af_"$af"_version_"$version"_seed_"$seed"_"
		python run_simulation.py $num_samples $num_tests $af $seed $version $output_file
	done
done


num_samples="2000"
num_tests="100"
af=".1"
seed="1"
version="none"
seed_array=("1" "2" "3" "4" "5")
for seed in "${seed_array[@]}"; do
	output_file=$simulation_results_dir"results_af_"$af"_version_"$version"_seed_"$seed"_"
	python run_simulation.py $num_samples $num_tests $af $seed $version $output_file
done
