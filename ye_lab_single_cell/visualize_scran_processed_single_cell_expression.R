library(ggplot2)
library(reshape2)
library(Matrix)
library(cowplot)
options(bitmapType='cairo')
options(warn=1)


figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=12), text = element_text(size=12),axis.text=element_text(size=10), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=10), legend.title = element_text(size=12)))
}




make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors <- function(categorical_variable, dim_reduction_1, dim_reduction_2, categorical_label, x_axis_label, y_axis_label) {
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

make_number_of_cells_per_individual_in_a_cell_type <- function(filtered_covariate_data, cell_type) {
  ct_indices <- filtered_covariate_data$cg_cov==cell_type
  total_cells = sum(ct_indices)
  ct_covariate_data <- filtered_covariate_data[ct_indices, ] 
  unique_indis <- as.character(unique(ct_covariate_data$ind_cov))
  cells_arr <- c()
  for (indi_index in 1:length(unique_indis)) {
    indi <- unique_indis[indi_index]
    num_cells <- sum(ct_covariate_data$ind_cov == indi)
    cells_arr <- c(cells_arr, num_cells)
  }
  ordering <- order(cells_arr)
  df <- data.frame(num_cells=cells_arr[ordering], indi=factor(unique_indis[ordering], levels=unique_indis[ordering]))

  p<-ggplot(data=df, aes(x=indi, y=num_cells)) +
    geom_bar(stat="identity", fill="steelblue")+
    figure_theme() +
    geom_hline(yintercept=50, linetype="dashed", color = "red") +
    geom_hline(yintercept=100, linetype="dashed", color = "green") +
    labs(x = "Individual", y = "Number of cells", title=paste0(cell_type, " (", total_cells, " total cells)")) +
      theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  return(p)
}

make_number_of_cells_per_individual_per_cell_type_bar_plots <- function(filtered_covariate_data) {
 b_bar_plot <- make_number_of_cells_per_individual_in_a_cell_type(filtered_covariate_data, "B")
 nk_bar_plot <- make_number_of_cells_per_individual_in_a_cell_type(filtered_covariate_data, "NK")
 pb_bar_plot <- make_number_of_cells_per_individual_in_a_cell_type(filtered_covariate_data, "PB")
 progen_bar_plot <- make_number_of_cells_per_individual_in_a_cell_type(filtered_covariate_data, "Progen")
 prolif_bar_plot <- make_number_of_cells_per_individual_in_a_cell_type(filtered_covariate_data, "Prolif")
 t4_bar_plot <- make_number_of_cells_per_individual_in_a_cell_type(filtered_covariate_data, "T4")
 t8_bar_plot <- make_number_of_cells_per_individual_in_a_cell_type(filtered_covariate_data, "T8")
 cdc_bar_plot <- make_number_of_cells_per_individual_in_a_cell_type(filtered_covariate_data, "cDC")
 cm_bar_plot <- make_number_of_cells_per_individual_in_a_cell_type(filtered_covariate_data, "cM")
 ncm_bar_plot <- make_number_of_cells_per_individual_in_a_cell_type(filtered_covariate_data, "ncM")
 pdc_bar_plot <- make_number_of_cells_per_individual_in_a_cell_type(filtered_covariate_data, "pDC")

  combined <- plot_grid(b_bar_plot, nk_bar_plot, pb_bar_plot, progen_bar_plot, prolif_bar_plot, t4_bar_plot, t8_bar_plot, cdc_bar_plot, cm_bar_plot, ncm_bar_plot, pdc_bar_plot, ncol = 2)
  return(combined)
}

make_number_of_cells_per_individual_bar_plot <- function(filtered_covariate_data) {
  unique_indis <- as.character(unique(filtered_covariate_data$ind_cov))
  cells_arr <- c()
  for (indi_index in 1:length(unique_indis)) {
    indi <- unique_indis[indi_index]
    num_cells <- sum(filtered_covariate_data$ind_cov == indi)
    cells_arr <- c(cells_arr, num_cells)
  }
  ordering <- order(cells_arr)
  df <- data.frame(num_cells=cells_arr[ordering], indi=factor(unique_indis[ordering], levels=unique_indis[ordering]))

  p<-ggplot(data=df, aes(x=indi, y=num_cells)) +
    geom_bar(stat="identity", fill="steelblue")+
    figure_theme() +
    labs(x = "Individual", y = "Number of cells") +
      theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  return(p) 
}

make_number_of_cells_per_cluster_in_a_cell_type <- function(pseudobulk_covariate_data, filtered_covariate_data, cluster_ids, cell_type) {
  cell_type_indices = pseudobulk_covariate_data$cg_cov_mode == cell_type
  unique_indis <- as.character(pseudobulk_covariate_data$pseudobulk_sample[cell_type_indices])
  
  cells_arr <- c()
  cluster_ids_string <- as.character(cluster_ids)
  for (indi_index in 1:length(unique_indis)) {
    indi <- unique_indis[indi_index]
    num_cells <- sum(cluster_ids_string == indi)
    cells_arr <- c(cells_arr, num_cells)
  }
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

make_number_of_cells_per_cluster_per_cell_type_bar_plot <- function(filtered_covariate_data, pseudobulk_covariate_data, cluster_ids) {
  b_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, filtered_covariate_data, cluster_ids, "B")
  nk_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, filtered_covariate_data, cluster_ids, "NK")
  pb_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, filtered_covariate_data, cluster_ids, "PB")
  progen_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, filtered_covariate_data, cluster_ids, "Progen")
  prolif_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, filtered_covariate_data, cluster_ids, "Prolif")
  t4_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, filtered_covariate_data, cluster_ids, "T4")
  t8_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, filtered_covariate_data, cluster_ids, "T8")
  cdc_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, filtered_covariate_data, cluster_ids, "cDC")
  cm_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, filtered_covariate_data, cluster_ids, "cM")
  ncm_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, filtered_covariate_data, cluster_ids, "ncM")
  pdc_bar_plot <- make_number_of_cells_per_cluster_in_a_cell_type(pseudobulk_covariate_data, filtered_covariate_data, cluster_ids, "pDC")


  combined <- plot_grid(b_bar_plot, nk_bar_plot, pb_bar_plot, progen_bar_plot, prolif_bar_plot, t4_bar_plot, t8_bar_plot, cdc_bar_plot, cm_bar_plot, ncm_bar_plot, pdc_bar_plot, ncol = 2)
  return(combined)
}




make_number_of_cells_per_cluster_bar_plot <- function(filtered_covariate_data, cluster_ids) {
  unique_indis <- as.character(unique(cluster_ids))
  cells_arr <- c()
  for (indi_index in 1:length(unique_indis)) {
    indi <- unique_indis[indi_index]
    num_cells <- sum(cluster_ids == indi)
    cells_arr <- c(cells_arr, num_cells)
  }
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



make_ct_proportion_per_individual_bar_plot_stratefied_by_sle_status <- function(filtered_covariate_data) {
  sle_cov_data = filtered_covariate_data[as.character(filtered_covariate_data$SLE_status)=="SLE",]
  control_cov_data = filtered_covariate_data[as.character(filtered_covariate_data$SLE_status)=="Healthy",]
  sle_plot <- make_ct_proportion_per_individual_bar_plot(sle_cov_data, "SLE")
  control_plot <- make_ct_proportion_per_individual_bar_plot(sle_cov_data, "Control")

  legendy <- get_legend(sle_plot + theme(legend.position="bottom"))
  combined <- plot_grid(sle_plot + theme(legend.position="none"), control_plot + theme(legend.position="none"), legendy, rel_heights=c(1,.2), ncol = 2)
  return(combined)

}



make_ct_proportion_per_individual_bar_plot <- function(filtered_covariate_data, title) {
  unique_indis <- as.character(unique(filtered_covariate_data$ind_cov))
  unique_cts <- as.character(unique(filtered_covariate_data$cg_cov))
  ct_ordering <- c("B", "NK", "PB", "Progen", "Prolif", "T4", "T8", "cDC", "cM", "ncM", "pDC")
  cells_arr <- c()
  for (indi_index in 1:length(unique_indis)) {
    indi <- unique_indis[indi_index]
    num_cells <- sum(filtered_covariate_data$ind_cov == indi)
    cells_arr <- c(cells_arr, num_cells)
  }
  ordering <- order(cells_arr)
  indi_ordering <- unique_indis[ordering]

  indi_arr <- c()
  ct_arr <- c()
  proportion_arr <- c()

  for (indi_index in 1:length(unique_indis)) {
    indi <- unique_indis[indi_index]
    total_cells <- sum(filtered_covariate_data$ind_cov == indi)
    for (ct_index in 1:length(unique_cts)) {
      ct <- unique_cts[ct_index]
      ct_cells <- sum((filtered_covariate_data$ind_cov==indi) & (filtered_covariate_data$cg_cov==ct))
      frac <- ct_cells/total_cells
      indi_arr <- c(indi_arr, indi)
      ct_arr <- c(ct_arr, ct)
      proportion_arr <- c(proportion_arr, frac)
    }
  }

  df <- data.frame(cell_type=factor(ct_arr, levels=ct_ordering), indi=factor(indi_arr, levels=indi_ordering), fraction=proportion_arr)

  plotter <- ggplot(data=df, aes(x=indi, y=fraction, fill=cell_type)) +
          geom_bar(stat="identity") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          figure_theme() +
          labs(x = "Individual", y = "Fraction", fill= "", title=title) +
          theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
  return(plotter)
}

make_batch_proportion_per_individual_bar_plot <- function(filtered_covariate_data) {
  unique_indis <- as.character(unique(filtered_covariate_data$ind_cov))
  unique_batches <- as.character(unique(filtered_covariate_data$batch_cov))
  cells_arr <- c()
  for (indi_index in 1:length(unique_indis)) {
    indi <- unique_indis[indi_index]
    #num_cells <- sum(filtered_covariate_data$ind_cov == indi)
    indi_batches = filtered_covariate_data$batch_cov[filtered_covariate_data$ind_cov == indi]
    mode = names(sort(summary(as.factor(indi_batches)), decreasing=T)[1])
    cells_arr <- c(cells_arr, mode)
  }
  ordering <- order(cells_arr)
  indi_ordering <- unique_indis[ordering]

  indi_arr <- c()
  ct_arr <- c()
  proportion_arr <- c()

  for (indi_index in 1:length(unique_indis)) {
    indi <- unique_indis[indi_index]
    total_cells <- sum(filtered_covariate_data$ind_cov == indi)
    for (ct_index in 1:length(unique_batches)) {
      ct <- unique_batches[ct_index]
      ct_cells <- sum((filtered_covariate_data$ind_cov==indi) & (filtered_covariate_data$batch_cov==ct))
      frac <- ct_cells/total_cells
      indi_arr <- c(indi_arr, indi)
      ct_arr <- c(ct_arr, ct)
      proportion_arr <- c(proportion_arr, frac)
    }
  }

  df <- data.frame(batch=factor(ct_arr), indi=factor(indi_arr, levels=indi_ordering), fraction=proportion_arr)

  plotter <- ggplot(data=df, aes(x=indi, y=fraction, fill=batch)) +
          geom_bar(stat="identity") +
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
          figure_theme() +
          labs(x = "Individual", y = "Fraction", fill= "Batch") +
          theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
  return(plotter)
}


######################################
# Make correlation heatmap correlating covariates with loadings
#######################################
make_covariate_loading_correlation_heatmap <- function(covariates, loadings) {
  loadings <- as.matrix(loadings)[,1:30]
  #print(summary(covariates))
  covariates$B = 1.0*(as.character(covariates$cg_cov) == "B")
  covariates$NK = 1.0*(as.character(covariates$cg_cov) == "NK")
  covariates$PB = 1.0*(as.character(covariates$cg_cov) == "PB")
  covariates$Progen = 1.0*(as.character(covariates$cg_cov) == "Progen")
  covariates$Prolif = 1.0*(as.character(covariates$cg_cov) == "Prolif")
  covariates$T4 = 1.0*(as.character(covariates$cg_cov) == "T4")
  covariates$T8 = 1.0*(as.character(covariates$cg_cov) == "T8")
  covariates$cDC = 1.0*(as.character(covariates$cg_cov) == "cDC")
  covariates$cM = 1.0*(as.character(covariates$cg_cov) == "cM")
  covariates$ncM = 1.0*(as.character(covariates$cg_cov) == "ncM")
  covariates$pDC = 1.0*(as.character(covariates$cg_cov) == "pDC")


  valid_covariates <- c(1, 2, 3, 5, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32)
  covariate_type <- c("cat", "cat", "cat", "cat", "num", "cat", "cat", "cat", "cat", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num")


  #print(length(valid_covariates))
  #print(length(covariate_type))

  num_cov = length(valid_covariates)

  cov_names <- colnames(covariates)[valid_covariates]
  num <- length(cov_names)
  # print(cov_names)



  covs <- covariates[,valid_covariates]
  print(head(covs))


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

######################################
# Make correlation heatmap correlating covariates with loadings
#######################################
make_pseudobulk_covariate_loading_correlation_heatmap <- function(covariates, loadings) {
  loadings <- as.matrix(loadings)[,1:50]


  valid_covariates <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)
  covariate_type <- c("cat", "num", "cat", "cat", "cat", "cat", "num", "cat", "cat", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num")

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

#########################
# Command line args
##########################
args <- commandArgs(TRUE)
processed_expression_dir <- args[1]  # Input dir
visualize_processed_expression_dir <- args[2]  # Output Dir
regress_out_batch <- args[3]  # Hyperparameter




##########################
# Load in Covariate data
##########################
# Load in Covariates
filtered_covariate_file <- paste0(processed_expression_dir, "scran_normalization_regress_batch_", regress_out_batch, "_cell_covariates.txt")
filtered_covariate_data <- read.table(filtered_covariate_file, header=TRUE, sep="\t")
saveRDS(filtered_covariate_data, "cov.rds")
#filtered_covariate_data <- readRDS("cov.rds")

# Load in PCS
pc_file <- paste0(processed_expression_dir, "scran_normalization_regress_batch_", regress_out_batch, "_sc_expression_pcs.txt")
pcs <- read.table(pc_file, header=FALSE, sep="\t")
saveRDS(pcs, "pcs.rds")
#pcs <- readRDS("pcs.rds")

# Load in PC PVE
pc_pve_file <- paste0(processed_expression_dir, "scran_normalization_regress_batch_", regress_out_batch, "_sc_expression_pcs_percent_variance_explained.txt")
pc_pve <- read.table(pc_pve_file, header=FALSE, sep="\t")

# Load in umap_loadings
umap_file <- paste0(processed_expression_dir, "scran_normalization_regress_batch_", regress_out_batch, "_sc_expression_umaps.txt")
umap_loadings <- read.table(umap_file, header=FALSE, sep="\t")
saveRDS(umap_loadings, "umap.rds")
#umap_loadings <- readRDS("umap.rds")

# Stem for output files
output_root <- paste0(visualize_processed_expression_dir, "regress_batch_", regress_out_batch, "_")

##########################
# Make bar plot showing number of cells per individual per cell type
##########################
num_cells_bar_plot <- make_number_of_cells_per_individual_per_cell_type_bar_plots(filtered_covariate_data)
output_file <- paste0(output_root, "number_of_cells_per_individual_per_cell_type_bar_plots.pdf")
ggsave(num_cells_bar_plot, file=output_file, width=7.2, height=10, units="in")

##########################
# Make bar plot showing number of cells per individual
##########################
num_cells_bar_plot <- make_number_of_cells_per_individual_bar_plot(filtered_covariate_data)
output_file <- paste0(output_root, "number_of_cells_per_individual_bar_plot.pdf")
ggsave(num_cells_bar_plot, file=output_file, width=7.2, height=5, units="in")

##########################
# Make stacked bar plot showing fraction of cell type per individual
##########################
ct_proportion_bar_plot <- make_ct_proportion_per_individual_bar_plot(filtered_covariate_data, "")
output_file <- paste0(output_root, "cell_type_proportions_per_individual_bar_plot.pdf")
ggsave(ct_proportion_bar_plot, file=output_file, width=7.2, height=5, units="in")

##########################
# Make stacked bar plot showing fraction of cell type per individual (stratefied by SLE status)
##########################
ct_proportion_bar_plot <- make_ct_proportion_per_individual_bar_plot_stratefied_by_sle_status(filtered_covariate_data)
output_file <- paste0(output_root, "cell_type_proportions_per_individual_stratefied_by_sle_status_bar_plot.pdf")
ggsave(ct_proportion_bar_plot, file=output_file, width=7.2, height=5, units="in")

##########################
# PCA-covariate heatmap
##########################
heatmap <- make_covariate_loading_correlation_heatmap(filtered_covariate_data, pcs)
output_file <- paste0(output_root, "covariate_pca_pve_heatmap.pdf")
ggsave(heatmap, file=output_file, width=7.2, height=10, units="in")


##########################
# Make PCA PVE line plot
##########################
num_pcs <- 50
output_file <- paste0(output_root, "pca_variance_explained_", num_pcs, "_pcs_line_plot.pdf")
ve_line_plot <- make_pc_variance_explained_line_plot(pc_pve[,1], num_pcs)
ggsave(ve_line_plot, file=output_file, width=7.2, height=5.0, units="in")



##########################
# Make PCA Plot colored by cell type
##########################
pca_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(filtered_covariate_data$cg_cov, pcs[,1], pcs[,2], "Cell Type", "PC1", "PC2")
output_file <- paste0(output_root, "pca_1_2_scatter_colored_by_cell_type.pdf")
ggsave(pca_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")


##########################
# Make PCA Plot colored by Lupus
##########################
pca_scatter_colored_by_lupus <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(filtered_covariate_data$SLE_status, pcs[,1], pcs[,2], "", "PC1", "PC2")
output_file <- paste0(output_root, "pca_1_2_scatter_colored_by_disease_cov.pdf")
ggsave(pca_scatter_colored_by_lupus, file=output_file, width=7.2, height=5, units="in")

##########################
# Make PCA Plot colored by Lupus
##########################
pca_scatter_colored_by_lupus <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(filtered_covariate_data$Status, pcs[,1], pcs[,2], "", "PC1", "PC2")
output_file <- paste0(output_root, "pca_1_2_scatter_colored_by_status.pdf")
ggsave(pca_scatter_colored_by_lupus, file=output_file, width=7.2, height=5, units="in")


##########################
# Make UMAP Plot colored by cell type
##########################
umap_scatter_colored_by_cell_type <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(filtered_covariate_data$cg_cov, umap_loadings[,1], umap_loadings[,2], "Cell Type", "umap1", "umap2")
output_file <- paste0(output_root, "umap_1_2_scatter_colored_by_cell_type.pdf")
ggsave(umap_scatter_colored_by_cell_type, file=output_file, width=7.2, height=5, units="in")

##########################
# Make UMAP Plot colored by status
##########################
umap_scatter_colored_by_batch <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(filtered_covariate_data$Status, umap_loadings[,1], umap_loadings[,2], "", "umap1", "umap2")
output_file <- paste0(output_root, "umap_1_2_scatter_colored_by_status.pdf")
ggsave(umap_scatter_colored_by_batch, file=output_file, width=7.2, height=5, units="in")

##########################
# Make UMAP Plot colored by Lupus
##########################
umap_scatter_colored_by_lupus <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(filtered_covariate_data$SLE_status, umap_loadings[,1], umap_loadings[,2], "", "umap1", "umap2")
output_file <- paste0(output_root, "umap_1_2_scatter_colored_by_disease_cov.pdf")
ggsave(umap_scatter_colored_by_lupus, file=output_file, width=7.2, height=5, units="in")

##########################
# Make UMAP Plot colored by Batch cov
##########################
umap_scatter_colored_by_batch <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(filtered_covariate_data$batch_cov, umap_loadings[,1], umap_loadings[,2], "", "umap1", "umap2")
output_file <- paste0(output_root, "umap_1_2_scatter_colored_by_batch_cov.pdf")
ggsave(umap_scatter_colored_by_batch, file=output_file, width=7.2, height=5, units="in")

##########################
# Make UMAP Plot colored by pop cov
##########################
umap_scatter_colored_by_pop <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(filtered_covariate_data$pop_cov, umap_loadings[,1], umap_loadings[,2], "", "umap1", "umap2")
output_file <- paste0(output_root, "umap_1_2_scatter_colored_by_pop_cov.pdf")
ggsave(umap_scatter_colored_by_pop, file=output_file, width=7.2, height=5, units="in")


##########################
# Make UMAP Plot colored by Ind cov
##########################
umap_scatter_colored_by_ind <- make_dimensionality_reduction_scatter_colored_by_categorical_variable_with_specified_cell_type_colors(filtered_covariate_data$ind_cov, umap_loadings[,1], umap_loadings[,2], "", "umap1", "umap2")
output_file <- paste0(output_root, "umap_1_2_scatter_colored_by_ind_cov.pdf")
ggsave(umap_scatter_colored_by_ind + theme(legend.position="none"), file=output_file, width=7.2, height=5, units="in")

















































