library(ggplot2)
library(reshape2)
library(Matrix)
library(cowplot)
options(bitmapType='cairo')
options(warn=1)


figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=12), text = element_text(size=12),axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=10), legend.title = element_text(size=12)))
}



make_dimensionality_reduction_scatter_colored_by_categorical_variable <- function(categorical_variable, dim_reduction_1, dim_reduction_2, categorical_label, x_axis_label, y_axis_label) {

	df <- data.frame(categorical_variable=factor(categorical_variable), dim_reduction_1=dim_reduction_1, dim_reduction_2=dim_reduction_2)
	

	scatter <- ggplot(df, aes(x=dim_reduction_1, y=dim_reduction_2, color=categorical_variable)) +
  				geom_point(size=.01) +
  				figure_theme() + 
  				labs(x=x_axis_label,y=y_axis_label, color=categorical_label) +
  				guides(colour = guide_legend(override.aes = list(size=3))) +
          theme(legend.position="none")
  	
  	
  	return(scatter)
}




make_pc_variance_explained_line_plot <- function(variance_explained, num_pcs) {
	variance_explained <- variance_explained[1:num_pcs]
	df <- data.frame(variance_explained = variance_explained, pc_num = 1:num_pcs)

	# PLOT AWAY
    line_plot <- ggplot(data=df, aes(x=pc_num, y=variance_explained)) +
                geom_line() +
                geom_point() +
                ylim(0,max(variance_explained) + .002) + 
                scale_x_continuous(breaks=seq(0,(num_pcs-1),5)) +
                labs(x = "PC number", y = "Variance Explained") + 
                figure_theme() 

    return(line_plot)
}


#########################
# Command line args
##########################
args <- commandArgs(TRUE)
processed_data_dir <- args[1]  # Input dir
visualize_processed_data_dir <- args[2]  # Output Dir



##########################
# Load in data
##########################
# Load in Covariates
cell_covariate_file <- paste0(processed_data_dir, "cell_covariates.txt")
covariate_data <- read.table(cell_covariate_file, sep="\t", header=TRUE, comment.char = "*")

# Expression pc loadings
pca_loadings_file <- paste0(processed_data_dir, "standardized_10_cap_normalized_expression_scanpy_hvg_all_genotyped_cells_pca_loadings.txt")
pca_loadings <- read.table(pca_loadings_file, header=FALSE, sep="\t")

# Load in PCA PVE
pca_pve_file <- paste0(processed_data_dir, "standardized_10_cap_normalized_expression_scanpy_hvg_all_genotyped_cells_pca_pve.txt")
pca_pve <- read.table(pca_pve_file, header=FALSE, sep="\t")


##########################
# Make PCA PVE line plot
##########################
num_pcs <- 50
output_file <- paste0(visualize_processed_data_dir, "pca_pve_line_plot_scanpy_hvg.pdf")
ve_line_plot <- make_pc_variance_explained_line_plot(pca_pve[,1], num_pcs)
ggsave(ve_line_plot, file=output_file, width=7.2, height=5.0, units="in")

##########################
# Make PCA Plot colored by differentiation
##########################
pc_scatter_colored_by_day <- make_dimensionality_reduction_scatter_colored_by_categorical_variable(covariate_data$day, pca_loadings[,1], pca_loadings[,2], "", "PC1", "PC2")
output_file <- paste0(visualize_processed_data_dir, "pc_scatter_colored_by_differentiation_day_scanpy_hvg.pdf")
ggsave(pc_scatter_colored_by_day, file=output_file, width=7.2, height=5, units="in")

pc_scatter_colored_by_plate_id <- make_dimensionality_reduction_scatter_colored_by_categorical_variable(factor(covariate_data$plate_id), pca_loadings[,1], pca_loadings[,2], "", "PC1", "PC2")
output_file <- paste0(visualize_processed_data_dir, "pc_scatter_colored_by_plate_id_scanpy_hvg.pdf")
ggsave(pc_scatter_colored_by_plate_id, file=output_file, width=7.2, height=5, units="in")

pc_scatter_colored_by_experiment <- make_dimensionality_reduction_scatter_colored_by_categorical_variable(factor(covariate_data$experiment), pca_loadings[,1], pca_loadings[,2], "", "PC1", "PC2")
output_file <- paste0(visualize_processed_data_dir, "pc_scatter_colored_by_experiment_scanpy_hvg.pdf")
ggsave(pc_scatter_colored_by_experiment, file=output_file, width=7.2, height=5, units="in")

pc_scatter_colored_by_donor <- make_dimensionality_reduction_scatter_colored_by_categorical_variable(factor(covariate_data$donor_long_id), pca_loadings[,1], pca_loadings[,2], "", "PC1", "PC2")
output_file <- paste0(visualize_processed_data_dir, "pc_scatter_colored_by_donor_scanpy_hvg.pdf")
ggsave(pc_scatter_colored_by_donor, file=output_file, width=7.2, height=5, units="in")





















##########################
# Load in data
##########################
# Load in Covariates
cell_covariate_file <- paste0(processed_data_dir, "cell_covariates.txt")
covariate_data <- read.table(cell_covariate_file, sep="\t", header=TRUE, comment.char = "*")

# Expression pc loadings
pca_loadings_file <- paste0(processed_data_dir, "standardized_10_cap_normalized_expression_scran_hvg_all_genotyped_cells_pca_loadings.txt")
pca_loadings <- read.table(pca_loadings_file, header=FALSE, sep="\t")

# Load in PCA PVE
pca_pve_file <- paste0(processed_data_dir, "standardized_10_cap_normalized_expression_scran_hvg_all_genotyped_cells_pca_pve.txt")
pca_pve <- read.table(pca_pve_file, header=FALSE, sep="\t")


##########################
# Make PCA PVE line plot
##########################
num_pcs <- 50
output_file <- paste0(visualize_processed_data_dir, "pca_pve_line_plot_scran_hvg.pdf")
ve_line_plot <- make_pc_variance_explained_line_plot(pca_pve[,1], num_pcs)
ggsave(ve_line_plot, file=output_file, width=7.2, height=5.0, units="in")

##########################
# Make PCA Plot colored by differentiation
##########################
pc_scatter_colored_by_day <- make_dimensionality_reduction_scatter_colored_by_categorical_variable(covariate_data$day, pca_loadings[,1], pca_loadings[,2], "", "PC1", "PC2")
output_file <- paste0(visualize_processed_data_dir, "pc_scatter_colored_by_differentiation_day_scran_hvg.pdf")
ggsave(pc_scatter_colored_by_day, file=output_file, width=7.2, height=5, units="in")

pc_scatter_colored_by_plate_id <- make_dimensionality_reduction_scatter_colored_by_categorical_variable(factor(covariate_data$plate_id), pca_loadings[,1], pca_loadings[,2], "", "PC1", "PC2")
output_file <- paste0(visualize_processed_data_dir, "pc_scatter_colored_by_plate_id_scran_hvg.pdf")
ggsave(pc_scatter_colored_by_plate_id, file=output_file, width=7.2, height=5, units="in")

pc_scatter_colored_by_experiment <- make_dimensionality_reduction_scatter_colored_by_categorical_variable(factor(covariate_data$experiment), pca_loadings[,1], pca_loadings[,2], "", "PC1", "PC2")
output_file <- paste0(visualize_processed_data_dir, "pc_scatter_colored_by_experiment_scran_hvg.pdf")
ggsave(pc_scatter_colored_by_experiment, file=output_file, width=7.2, height=5, units="in")

pc_scatter_colored_by_donor <- make_dimensionality_reduction_scatter_colored_by_categorical_variable(factor(covariate_data$donor_long_id), pca_loadings[,1], pca_loadings[,2], "", "PC1", "PC2")
output_file <- paste0(visualize_processed_data_dir, "pc_scatter_colored_by_donor_scran_hvg.pdf")
ggsave(pc_scatter_colored_by_donor, file=output_file, width=7.2, height=5, units="in")





