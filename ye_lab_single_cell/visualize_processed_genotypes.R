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
	# Put everthing in compact data frame
	unique_categories = as.character(unique(categorical_variable))

	df <- data.frame(categorical_variable=factor(categorical_variable, levels=unique_categories), dim_reduction_1=dim_reduction_1, dim_reduction_2=dim_reduction_2)
	
	scatter <- ggplot(df, aes(x=dim_reduction_1, y=dim_reduction_2, color=categorical_variable)) +
  				geom_point(size=1) +
  				figure_theme() + 
  				labs(x=x_axis_label,y=y_axis_label, color=categorical_label) +
  				guides(colour = guide_legend(override.aes = list(size=3)))
  	return(scatter)
}








#########################
# Command line args
##########################
args <- commandArgs(TRUE)
processed_genotype_dir <- args[1]  # Input dir
visualize_processed_genotype_dir <- args[2]  # Output Dir




# Genotype pc results file
genotype_pcs_file <- paste0(processed_genotype_dir, "individual_level_genotype_pcs.txt")
genotype_pcs <- read.table(genotype_pcs_file, header=TRUE, sep="\t")

##########################
# Make PCA Plot colored by ancestry
##########################
pca_scatter_colored_by_ancestry <- make_dimensionality_reduction_scatter_colored_by_categorical_variable(genotype_pcs$ancestry, genotype_pcs$geno_pc1, genotype_pcs$geno_pc2, "", "Genotype PC1", "Genotype PC2")
output_file <- paste0(visualize_processed_genotype_dir, "geno_pca_1_2_scatter_colored_by_ancestry.pdf")
ggsave(pca_scatter_colored_by_ancestry, file=output_file, width=7.2, height=5, units="in")