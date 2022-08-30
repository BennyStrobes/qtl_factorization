




#############################
# Output data
#############################

# Simulation output root
sim_output_root="/scratch16/abattle4/bstrober/qtl_factorization/simulation/"

# Data containing simulated data
sim_data_dir=$sim_output_root"simulated_data/"

# Directory containing eqtl factorization results run on simulated data
eqtl_results_dir=$sim_output_root"surge_results/"

# Directory containing visualizations of the results
viz_dir=$sim_output_root"visualize_simulation/"




sh run_simulation_experiments.sh $sim_data_dir $eqtl_results_dir