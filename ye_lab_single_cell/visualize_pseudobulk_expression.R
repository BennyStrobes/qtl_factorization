library(ggplot2)
library(reshape2)
library(Matrix)
library(cowplot)
options(bitmapType='cairo')
options(warn=1)


figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=12), text = element_text(size=12),axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=10), legend.title = element_text(size=12)))
}

make_knn_cell_type_stacked_bar_chart <- function(knn_ct_summary_df) {
  cell_type_vec <- c()
  neighbor_cell_type_vec <- c()
  fraction_vec <- c()
  num_cell_types <- dim(knn_ct_summary_df)[1]
  cell_types <- as.character(knn_ct_summary_df$cell_type)

  for (primary_ct_num in 1:num_cell_types) {
    primary_ct <- cell_types[primary_ct_num]
    for (knn_ct_num in 1:num_cell_types) {
      knn_ct <- cell_types[knn_ct_num]
      frac <- knn_ct_summary_df[primary_ct_num, (knn_ct_num + 1)]
      cell_type_vec <- c(cell_type_vec, primary_ct)
      neighbor_cell_type_vec <- c(neighbor_cell_type_vec, knn_ct)
      fraction_vec <- c(fraction_vec, frac)
    }
  }

  df <- data.frame(cell_type=factor(cell_type_vec, levels=cell_types), knn_cell_type=factor(neighbor_cell_type_vec, levels=cell_types), fraction=fraction_vec)

  plotter <- ggplot(data=df, aes(x=cell_type, y=fraction, fill=knn_cell_type)) +
          geom_bar(stat="identity") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          figure_theme()
  return(plotter)
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

######################################
# Make correlation heatmap correlating covariates with loadings
#######################################
make_pseudobulk_covariate_loading_correlation_heatmap <- function(covariates, loadings) {
  loadings <- as.matrix(loadings)[,1:50]


  valid_covariates <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22)
  covariate_type <- c("cat", "num", "cat", "cat", "cat", "cat", "num", "cat", "cat", "cat", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num")

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




make_number_of_cells_per_cluster_bar_plot <- function(pseudobulk_covariate_data) {
  cells_arr <- pseudobulk_covariate_data$num_cells
  unique_indis <- pseudobulk_covariate_data$pseudobulk_sample
  ordering <- order(cells_arr)
  df <- data.frame(num_cells=cells_arr[ordering], indi=factor(unique_indis[ordering], levels=unique_indis[ordering]))

  p<-ggplot(data=df, aes(x=indi, y=num_cells)) +
    geom_bar(stat="identity", fill="steelblue")+
    figure_theme() +
    labs(x = "Clusters", y = "Number of cells") +
      theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  return(p) 
}

make_number_of_clusters_per_individual_bar_plot <- function(pseudobulk_covariate_data) {
  individuals <- as.character(pseudobulk_covariate_data$ind_cov)
  unique_indis <- as.character(unique(individuals))
  cells_arr <- c()
  for (indi_index in 1:length(unique_indis)) {
    indi <- unique_indis[indi_index]
    num_cells <- sum(individuals == indi)
    cells_arr <- c(cells_arr, num_cells)
  }
  ordering <- order(cells_arr)
  df <- data.frame(num_cells=cells_arr[ordering], indi=factor(unique_indis[ordering], levels=unique_indis[ordering]))

  p<-ggplot(data=df, aes(x=indi, y=num_cells)) +
    geom_bar(stat="identity", fill="steelblue")+
    figure_theme() +
    labs(x = "Individual", y = "Number of clusters") +
      theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  return(p) 
}

make_number_of_cells_per_cluster_in_a_cell_type <- function(pseudobulk_covariate_data, cell_type) {
  cell_type_indices = pseudobulk_covariate_data$cg_cov_mode == cell_type
  unique_indis <- as.character(pseudobulk_covariate_data$pseudobulk_sample[cell_type_indices])
 

  cells_arr <- pseudobulk_covariate_data$num_cells[cell_type_indices]
  ordering <- order(cells_arr)
  df <- data.frame(num_cells=cells_arr[ordering], indi=factor(unique_indis[ordering], levels=unique_indis[ordering]))

  p<-ggplot(data=df, aes(x=indi, y=num_cells)) +
    geom_bar(stat="identity", fill="steelblue")+
    figure_theme() +
    labs(x = "Clusters", y = "Number of cells", title=paste0(cell_type)) +
      theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
    geom_hline(yintercept=50, linetype="dashed", color = "red") +
    geom_hline(yintercept=100, linetype="dashed", color = "green")
  return(p) 


}


make_number_of_cells_per_cluster_per_cell_type_bar_plot <- function(pseudobulk_covariate_data) {
  b_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, "B")
  nk_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, "NK")
  pb_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, "PB")
  progen_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, "Progen")
  prolif_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, "Prolif")
  t4_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, "T4")
  t8_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, "T8")
  cdc_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, "cDC")
  cm_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, "cM")
  ncm_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, "ncM")
  pdc_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, "pDC")


  combined <- plot_grid(b_bar_plot, nk_bar_plot, pb_bar_plot, progen_bar_plot, prolif_bar_plot, t4_bar_plot, t8_bar_plot, cdc_bar_plot, cm_bar_plot, ncm_bar_plot, pdc_bar_plot, ncol = 2)
  return(combined)
}

make_dimensionality_reduction_scatter_colored_by_categorical_variable <- function(categorical_variable, dim_reduction_1, dim_reduction_2, categorical_label, x_axis_label, y_axis_label) {
	# Put everthing in compact data frame
	unique_categories = as.character(unique(categorical_variable))

	df <- data.frame(categorical_variable=factor(categorical_variable, levels=unique_categories), dim_reduction_1=dim_reduction_1, dim_reduction_2=dim_reduction_2)
	

	scatter <- ggplot(df, aes(x=dim_reduction_1, y=dim_reduction_2, color=categorical_variable)) +
  				geom_point(size=.01) +
  				figure_theme() + 
  				labs(x=x_axis_label,y=y_axis_label, color=categorical_label) +
  				guides(colour = guide_legend(override.aes = list(size=3)))
  	
  	
  	return(scatter)
}






#########################
# Command line args
##########################
args <- commandArgs(TRUE)
processed_pseudobulk_expression_dir <- args[1]  # input dir
cluster_resolution <- args[2]  # Hyperparameter
visualize_processed_pseudobulk_expression_dir <- args[3]  # output dir
regress_out_batch <- args[4]
gene_level_normalization <- args[5]


############################
# Load in data
############################
# Get cluster neighbor ct summary file
neighbor_ct_summary_file <- paste0(processed_pseudobulk_expression_dir, "scran_normalization_regress_batch_", regress_out_batch, "_individual_clustering_leiden_resolution_", cluster_resolution, "_cell_type_summary.txt")
neighbor_ct_summary_df <- read.table(neighbor_ct_summary_file, header=TRUE, sep="\t")

# Load in pseudobulk covariate data
pseudobulk_covariate_file <- paste0(processed_pseudobulk_expression_dir, "pseudobulk_scran_normalization_regress_batch_", regress_out_batch, "_individual_clustering_leiden_resolution_", cluster_resolution, "_sample_covariates.txt")
pseudobulk_covariate_data <- read.table(pseudobulk_covariate_file, header=TRUE, sep="\t")

# Load in pseudobulk expression_pcs
pseudobulk_pcs_file <- paste0(processed_pseudobulk_expression_dir, "pseudobulk_scran_normalization_regress_batch_", regress_out_batch, "_individual_clustering_leiden_resolution_", cluster_resolution, "_none_sample_norm_", gene_level_normalization, "_gene_norm_pca_scores.txt")
pseudobulk_pcs <- read.table(pseudobulk_pcs_file, header=FALSE, sep="\t")

# Load in pseudobulk PC PVE
pseudobulk_pc_pve_file <- paste0(processed_pseudobulk_expression_dir, "pseudobulk_scran_normalization_regress_batch_", regress_out_batch, "_individual_clustering_leiden_resolution_", cluster_resolution, "_none_sample_norm_", gene_level_normalization, "_gene_norm_pca_pve.txt")
pseudobulk_pc_pve <- read.table(pseudobulk_pc_pve_file, header=FALSE, sep="\t")


############################
# Output root
############################
output_root <- paste0(visualize_processed_pseudobulk_expression_dir, "pseudobulk_scran_normalization_regress_batch_", regress_out_batch, "_resolution_", cluster_resolution, "_gene_level_normalization_", gene_level_normalization, "_")

##########################
# Clustering neighboring cell type summary
##########################
cluster_ct_stacked_bar <- make_knn_cell_type_stacked_bar_chart(neighbor_ct_summary_df)
output_file <- paste0(output_root, "cluster_neighbor_cell_type_stacked_bar_chart.pdf")
ggsave(cluster_ct_stacked_bar, file=output_file, width=7.2, height=5, units="in")



##########################
# PCA-covariate heatmap for pseudobulk data
##########################
heatmap <- make_pseudobulk_covariate_loading_correlation_heatmap(pseudobulk_covariate_data, pseudobulk_pcs)
output_file <- paste0(output_root, "pseudobulk_covariate_pca_pve_heatmap.pdf")
ggsave(heatmap, file=output_file, width=7.2, height=10, units="in")


##########################
# Make PCA PVE line plot
##########################
num_pcs <- 50
output_file <- paste0(output_root, "pseudobulk_pca_variance_explained_", num_pcs, "_pcs_line_plot.pdf")
ve_line_plot <- make_pc_variance_explained_line_plot(pseudobulk_pc_pve[,1], num_pcs)
ggsave(ve_line_plot, file=output_file, width=7.2, height=5.0, units="in")



##########################
# Make bar plot showing number of cells per cluster
##########################
num_cells_bar_plot <- make_number_of_cells_per_cluster_bar_plot(pseudobulk_covariate_data)
output_file <- paste0(output_root, "number_of_cells_per_pseudobulk_cluster_bar_plot.pdf")
ggsave(num_cells_bar_plot, file=output_file, width=7.2, height=5, units="in")


##########################
# Make Number of clusters per individual
##########################
num_cells_bar_plot <- make_number_of_clusters_per_individual_bar_plot(pseudobulk_covariate_data)
output_file <- paste0(output_root, "number_of_clusters_per_individual_bar_plot.pdf")
ggsave(num_cells_bar_plot, file=output_file, width=7.2, height=5, units="in")

##########################
# Make bar plot showing number of cells per cluster per cell type
##########################
num_cells_bar_plot <- make_number_of_cells_per_cluster_per_cell_type_bar_plot(pseudobulk_covariate_data)
output_file <- paste0(output_root, "number_of_cells_per_pseudobulk_cluster_per_cell_type_bar_plot.pdf")
ggsave(num_cells_bar_plot, file=output_file, width=7.2, height=10, units="in")


##########################
# Make PCA Plot colored by cell type
##########################
pca_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable(pseudobulk_covariate_data$cg_cov_mode, pseudobulk_pcs[,1], pseudobulk_pcs[,2], "Cell Type", "PC1", "PC2")
output_file <- paste0(output_root, "pca_1_2_scatter_colored_by_cell_type.pdf")
ggsave(pca_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")

##########################
# Make PCA Plot colored by batch
##########################
pca_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable(pseudobulk_covariate_data$batch_cov, pseudobulk_pcs[,1], pseudobulk_pcs[,2], "Batch", "PC1", "PC2")
output_file <- paste0(output_root, "pca_1_2_scatter_colored_by_batch.pdf")
ggsave(pca_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")


##########################
# Make PCA Plot colored by Indi
##########################
pca_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable(pseudobulk_covariate_data$ind_cov, pseudobulk_pcs[,1], pseudobulk_pcs[,2], "Ind", "PC1", "PC2")
output_file <- paste0(output_root, "pca_1_2_scatter_colored_by_ind.pdf")
ggsave(pca_scatter_colored_by_cell_type + theme(legend.position="none"), file=output_file, width=7.2, height=5, units="in")

