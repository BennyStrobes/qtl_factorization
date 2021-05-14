args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(cowplot)
library(umap)
library(ggplot2)
library(RColorBrewer)
options(bitmapType = 'cairo', device = 'pdf')



gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}




make_loading_boxplot_plot_by_categorical_covariate <- function(covariates, loadings, covariate_name) {
	loading_vec <- c()
	covariate_vec <- c()
	factor_vec <- c()

	num_factors <- dim(loadings)[2]

	for (factor_number in 1:num_factors) {
		loading_vec <- c(loading_vec, loadings[,factor_number])
		covariate_vec <- c(covariate_vec, as.character(covariates))
		factor_vec <- c(factor_vec, rep(as.character(factor_number), length(covariates)))
	}


	df <- data.frame(loading=loading_vec, covariate=factor(covariate_vec), latent_factor=factor(factor_vec, levels=as.character(1:num_factors)))



	boxplot <- ggplot(df, aes(x=latent_factor, y=loading, fill=covariate)) + geom_boxplot(outlier.size = .00001) +
				gtex_v8_figure_theme() + 
	        	labs(x="Latent factor", y = "Sample loading", fill=covariate_name) +
	        	theme(legend.position="bottom") +
	        	guides(colour = guide_legend(override.aes = list(size=2))) +
	           	guides(colour=guide_legend(nrow=4,byrow=TRUE, override.aes = list(size=2)))
	   

	return(boxplot)
}



make_umap_loading_scatter_plot_colored_by_categorical_variable <- function(covariates, umap_loadings, covariate_name) {
	df <- data.frame(loading_1=umap_loadings[,1], loading_2=umap_loadings[,2], covariate=factor(covariates))
	plotter <- ggplot(df) + 
	           geom_point(aes(x=loading_1, y=loading_2, color=covariate), size=.001) +
	           gtex_v8_figure_theme() + 
	           guides(colour = guide_legend(override.aes = list(size=2))) +
	           labs(x="UMAP 1", y = "UMAP 2", color=covariate_name) + 
	           guides(colour=guide_legend(nrow=3,byrow=TRUE, override.aes = list(size=2))) +
	           theme(legend.position="bottom") + 
	           theme(legend.text = element_text(size=8), legend.title = element_text(size=8))
	return(plotter)
}


make_umap_loading_scatter_plot_colored_by_real_valued_variable <- function(covariates, umap_loadings, covariate_name) {
	df <- data.frame(loading_1=umap_loadings[,1], loading_2=umap_loadings[,2], covariate=covariates)
	plotter <- ggplot(df) + 
	           geom_point(aes(x=loading_1, y=loading_2, color=covariate), size=.01) +
	           gtex_v8_figure_theme() + 
	           labs(x="UMAP 1", y = "UMAP 2", color=covariate_name) + 
	           scale_color_gradient(low="pink",high="blue") +
	           theme(legend.position="bottom") + 
	           theme(legend.text = element_text(size=8), legend.title = element_text(size=8))
	return(plotter)
}


######################################
# Make correlation heatmap correlating covariates with loadings
#######################################
make_pseudobulk_covariate_loading_correlation_heatmap <- function(covariates, loadings) {
  loadings <- as.matrix(loadings)


  valid_covariates <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)
  covariate_type <- c("cat", "num", "cat", "cat", "cat", "cat", "num", "cat", "cat", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num")

  #print(summary(covariates))

  #print(length(valid_covariates))
  #print(length(covariate_type))

  num_cov = length(valid_covariates)

  cov_names <- colnames(covariates)[valid_covariates]
  num <- length(cov_names)
  # print(cov_names)



  covs <- covariates[,valid_covariates]


  # Initialize PVE heatmap
  factor_colnames <- paste0("Factor", 1:(dim(loadings)[2]))
  factor_rownames <- colnames(covs)
  pve_map <- matrix(0, dim(covs)[2], dim(loadings)[2])
  colnames(pve_map) <- factor_colnames
  rownames(pve_map) <- colnames(covs)


    # Loop through each PC, COV Pair and take correlation
    num_pcs <- dim(loadings)[2]
    num_covs <- dim(covs)[2]
    for (num_pc in 1:num_pcs) {
        for (num_cov in 1:num_covs) {
            pc_vec <- loadings[,num_pc]
            cov_vec <- covs[,num_cov]
            #print(paste0(num_pc, " - ", num_cov))
            if (covariate_type[num_cov] == "cat") {
            #print(cov_vec[1:10])
              lin_model <- lm(pc_vec ~ factor(cov_vec))
          } else {
            lin_model <- lm(pc_vec ~ cov_vec)
          }
            pve_map[num_cov, num_pc] <- summary(lin_model)$adj.r.squared
        }
    }
    
    ord <- hclust( dist(scale(pve_map), method = "euclidean"), method = "ward.D" )$order

    melted_mat <- melt(pve_map)
    colnames(melted_mat) <- c("Covariate", "Loading","PVE")

    melted_mat$Covariate = factor(melted_mat$Covariate, levels=rownames(pve_map)[ord])
    melted_mat$Loading = factor(melted_mat$Loading, levels=factor_colnames)
   #  Use factors to represent covariate and pc name
    # melted_mat$Covariate 
    # melted_mat$Covariate <- factor(melted_mat$Covariate, levels = rownames(pve_map)[ord])
    #melted_mat$PC <- substr(as.character(melted_mat$PC),3,5)
    #melted_mat$PC <- factor(melted_mat$PC, levels=paste0("", 1:(length(unique(melted_mat$PC)))))

    
    #levels(melted_mat$PC) = paste0("PC", 1:(length(levels(melted_mat$PC))))
    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=Covariate, y=Loading)) + geom_tile(aes(fill=PVE)) + scale_fill_gradient2(midpoint=-.05, guide="colorbar")
    heatmap <- heatmap + labs(y="",fill="VE")
    heatmap <- heatmap + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
    # Save File
    return(heatmap)
}








############################
# Command line args
############################
processed_data_dir <- args[1]
eqtl_results_dir <- args[2]
visualization_dir <- args[3]

############################
# Load in files
############################
sample_covariate_file <- paste0(processed_data_dir, "cluster_tmm_ign_pseudobulk_leiden_no_cap_2.5_sample_covariates.txt")

############################
# Model Specification
############################
model_stem <- paste0("eqtl_factorization_results_lf_interaction_egenes_no_cap_2.5_eqtl_factorization_vi_results_k_init_20_lambda_v_1_seed_1_var_param_1e-3_temper_")
eqtl_factorization_loading_file <- paste0(eqtl_results_dir, model_stem, "U_S.txt")

print(sample_covariate_file)
print(eqtl_factorization_loading_file)
# Load in data
covariates <- read.table(sample_covariate_file, header=TRUE, sep="\t")
loadings <- read.table(eqtl_factorization_loading_file, header=FALSE)

######################################
# Make loading boxplot colored by Ancestry
#######################################
output_file <- paste0(visualization_dir, "loading_boxplot_colored_by_ancestry.pdf")
boxplot <- make_loading_boxplot_plot_by_categorical_covariate(covariates$pop_cov, loadings, "Known ancestry")
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")

######################################
# Make loading boxplot colored by Ancestry
#######################################
output_file <- paste0(visualization_dir, "loading_boxplot_colored_by_cell_type.pdf")
boxplot <- make_loading_boxplot_plot_by_categorical_covariate(covariates$cg_cov_mode, loadings, "Cell Type")
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")

######################################
# Make loading boxplot colored by Sex
#######################################
output_file <- paste0(visualization_dir, "loading_boxplot_colored_by_sex.pdf")
boxplot <- make_loading_boxplot_plot_by_categorical_covariate(covariates$Sex, loadings, "Sex")
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")




##########################
# PCA-covariate heatmap for pseudobulk data
##########################
heatmap <- make_pseudobulk_covariate_loading_correlation_heatmap(covariates, loadings)
output_file <- paste0(visualization_dir, "eqtl_factorization_covariate_pca_pve_heatmap.pdf")
ggsave(heatmap, file=output_file, width=7.2, height=10, units="in")



if (FALSE) {

print('UMAP START')
umap_loadings = umap(loadings)$layout
#saveRDS( umap_loadings, "umap_loadings_k_90_median_gaussian_kernel.rds")
print('UMAP DONE')


######################################
# Visualize UMAP scatter plot colored by individual
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_individual.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$ind_cov, umap_loadings, "Known individual")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by cell type
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_cell_type.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$cg_cov_mode, umap_loadings, "cell type")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by number of cells
#######################################
output_file <- paste0(eqtl_visualization_dir, "umap_loading_scatter_colored_by_num_cells.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$num_cells, umap_loadings, "Number of cells")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")
}
