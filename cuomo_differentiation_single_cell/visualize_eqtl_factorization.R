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

make_loading_boxplot_for_one_factor_by_categorical_covariate <- function(covariate, loadings, factor_number, covariate_name) {
  df <- data.frame(loading=loadings, covariate=factor(covariate))

  boxplot <- ggplot(df, aes(x=covariate, y=loading, fill=covariate)) + geom_boxplot(outlier.size = .00001) +
        gtex_v8_figure_theme() + 
            labs(x="", y = paste0("Sample loading (", factor_number,")"), fill="") +
            theme(legend.position="none") +
            guides(colour = guide_legend(override.aes = list(size=2))) +
            theme(axis.text.x=element_blank()) + 
              guides(colour=guide_legend(nrow=4,byrow=TRUE, override.aes = list(size=2))) 

}



make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate <- function(covariate, loadings, covariate_name) {
  loading_vec <- c()
  covariate_vec <- c()
  num_factors <- dim(loadings)[2]

  plot_arr <- list()

  for (factor_num in 1:num_factors) {
    factor_boxplot <- make_loading_boxplot_for_one_factor_by_categorical_covariate(covariate, loadings[,factor_num], factor_num, covariate_name)
    plot_arr[[factor_num]] <- factor_boxplot
  }

  merged = plot_grid(plotlist=plot_arr, ncol=1)


  return(merged)
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

make_scatter_plot_colored_by_categorical_variable <- function(x_var, y_var, color_var, x_axis_label, y_axis_label, color_axis_label, title_label) {
  df <- data.frame(loading_1=x_var, loading_2=y_var, covariate=factor(color_var))
  plotter <- ggplot(df) + 
             geom_point(aes(x=loading_1, y=loading_2, color=covariate), size=.001) +
             gtex_v8_figure_theme() + 
             guides(colour = guide_legend(override.aes = list(size=2))) +
             labs(x=x_axis_label, y = y_axis_label, color=color_axis_label, title=title_label) + 
             guides(colour=guide_legend(nrow=3,byrow=TRUE, override.aes = list(size=2))) +
             theme(legend.position="bottom") + 
             theme(legend.text = element_text(size=8), legend.title = element_text(size=8))
  return(plotter)
}

make_loading_vs_loading_scatter_colored_by_categorical_covariate <- function(loadings1, loadings2, covariates, loading1_name, loading2_name, covariate_name) {
  df <- data.frame(loading_1=loadings1, loading_2=loadings2, covariate=factor(covariates))
  plotter <- ggplot(df) + 
             geom_point(aes(x=loading_1, y=loading_2, color=covariate), size=.01) +
             gtex_v8_figure_theme() + 
             guides(colour = guide_legend(override.aes = list(size=2))) +
             labs(x=loading1_name, y = loading2_name, color=covariate_name) + 
             guides(colour=guide_legend(nrow=3,byrow=TRUE, override.aes = list(size=2))) +
             theme(legend.position="bottom") + 
             theme(legend.text = element_text(size=8), legend.title = element_text(size=8))
  return(plotter)
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


make_loading_vs_loading_scatter_colored_by_real_valued_covariate <- function(loadings1, loadings2, covariates, loading1_name, loading2_name, covariate_name) {
  df <- data.frame(loading_1=loadings1, loading_2=loadings2, covariate=(covariates))
  plotter <- ggplot(df) + 
             geom_point(aes(x=loading_1, y=loading_2, color=covariate), size=.01) +
             gtex_v8_figure_theme() + 
             labs(x=loading1_name, y = loading2_name, color=covariate_name) + 
             scale_color_gradient(low="pink",high="blue") +
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


  valid_covariates <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)
  covariate_type <- c("cat", "num", "cat", "cat", "cat", "cat", "num", "cat", "cat", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num")

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


make_factor_distribution_histograms <- function(factors) {
  factor_min <- floor(min(factors)) 
  factor_max <- ceiling(max(factors))

  num_factors <- dim(factors)[1]

  plot_arr <- list()
  factor_vals_vec <- c()
  factor_names_vec <- c()
  factor_num <- 1
  num_tests <- dim(factors)[2]
  print(num_factors)
  ordered_factors <- c()
  for (factor_num in 1:num_factors) {
    factor_vals_vec <- c(factor_vals_vec, as.numeric(factors[factor_num,]))
    factor_names_vec <- c(factor_names_vec, rep(paste0("factor", factor_num), num_tests))
    ordered_factors <- c(ordered_factors, paste0("factor", factor_num))
  }
  df <- data.frame(factor_values=factor_vals_vec, factor_names=factor(factor_names_vec, levels=ordered_factors))
  p <- ggplot(df, aes(x=factor_values))+
      geom_histogram(color="darkblue", fill="lightblue") + 
      facet_grid(factor_names ~ .)
      gtex_v8_figure_theme() 

  return(p)
}

make_pc_variance_explained_line_plot <- function(variance_explained) {
  num_pcs <- length(variance_explained)
  variance_explained <- variance_explained[1:num_pcs]
  df <- data.frame(variance_explained = variance_explained, pc_num = 1:num_pcs)

  # PLOT AWAY
    line_plot <- ggplot(data=df, aes(x=pc_num, y=variance_explained)) +
                geom_line() +
                geom_point() +
                ylim(0,max(variance_explained) + (max(variance_explained)*.1)) + 
                scale_x_continuous(breaks=seq(0,(num_pcs),1)) +
                labs(x = "Factor number", y = "Fraction of QTL variance explained") + 
                gtex_v8_figure_theme() 

    return(line_plot)
}

make_histogram_of_loadings_for_cell_type_stratefied_by_sle_status <- function(loadings, sle_status, cell_type, loading_num) {
  df <- data.frame(loadings=loadings, sle_status=factor(sle_status))

  # Wilcoxon rank sum test to investigate difference in distributions
  healthy_indices = as.character(sle_status) == "Healthy"
  sle_indices = as.character(sle_status) == "SLE"
  pvalue = wilcox.test(loadings[healthy_indices], loadings[sle_indices])$p.value

  p <- ggplot(data=df, aes(x=loadings, fill=sle_status)) +
    #geom_histogram(aes(y=..ncount..),color="#e9ecef", alpha=0.6, position = 'identity') +
    geom_density( color="#e9ecef", alpha=0.6, position = 'identity') +
    scale_fill_manual(values=c("#69b3a2", "#404080")) +
    gtex_v8_figure_theme() +
    labs(x=paste0("Factor ", loading_num, " loadings"), fill="", title=paste0(cell_type , " / wilcoxon pvalue: ", signif(pvalue, digits = 4)))
  return(p)

}

make_histogram_of_loadings_for_each_cell_type_stratefied_by_sle_status <- function(loadings, cg_cov, sle_status, loading_num) {
  cell_types <- sort(as.character(unique(cg_cov)))

  plot_arr <- list()
  for (cell_type_iter in 1:length(cell_types)) {
    cell_type <- cell_types[cell_type_iter]
    cell_type_indices = as.character(cg_cov) == cell_type
    histy <- make_histogram_of_loadings_for_cell_type_stratefied_by_sle_status(loadings[cell_type_indices], sle_status[cell_type_indices], cell_type, loading_num)
    plot_arr[[cell_type_iter]] <- histy + theme(legend.position="none")
  }
  legend <- get_legend(histy)
  plot_arr[[(cell_type_iter+1)]] = legend

  merged = plot_grid(plotlist=plot_arr, ncol=2)

  return(merged)
}

make_covariate_loading_correlation_heatmap <- function(covariates, loadings) {
  loadings <- as.matrix(loadings)


  valid_covariates <- c(7, 8, 10, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 31, 32, 33, 38, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 69,71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81,83, 85, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103)
  covariate_type <- c("cat", "cat", "cat", "real", "real", "real", "real", "real", "real", "real", "real", "real", "real", "real", "real", "real", "real", "real","real","real","real","real","real","real","real","real","real","real","real","real","real","real","real", "cat", "real", "real", "real","real","real","real","real","real","real","real","real","real","cat", "cat", "real","real","real","real","real","real","real","real","real","real","real","real","real","real","real","real","real")


  #print(summary(covariates))

  #print(length(valid_covariates))
  #print(length(covariate_type))

  num_cov = length(valid_covariates)

  cov_names <- colnames(covariates)[valid_covariates]

  #print(cov_names)
  #print(length(valid_covariates))
  #print(length(covariate_type))
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

make_loading_expression_pc_scatter_for_each_factor <- function(expression_pc, loading, factor_number, x_axis_label) {
  df <- data.frame(expression_pc=expression_pc, loading=loading)
  print(factor_number)
  print(cor(expression_pc, loading)*cor(expression_pc, loading))

  plotter <- ggplot(df, aes(x=expression_pc, y=loading)) + 
             geom_point(size=.001, alpha=.35) +
             gtex_v8_figure_theme() + 
             labs(x=x_axis_label, y = paste0("Loading ", factor_number)) + 
             theme(legend.text = element_text(size=8), legend.title = element_text(size=8)) +
             geom_smooth()
  return(plotter)

}


loading_expression_pc1_scatter_with_row_for_every_factor <- function(expression_pc, loadings, x_axis_label) {
  num_factors <- dim(loadings)[2]

  plot_arr <- list()

  for (factor_num in 1:num_factors) {
    factor_boxplot <- make_loading_expression_pc_scatter_for_each_factor(expression_pc, loadings[, factor_num], factor_num, x_axis_label)
    plot_arr[[factor_num]] <- factor_boxplot
  }

  merged = plot_grid(plotlist=plot_arr, ncol=1)


  return(merged)
}

############################
# Command line args
############################
processed_data_dir <- args[1]
eqtl_results_dir <- args[2]
visualization_dir <- args[3]
model_stem <- args[4]
output_stem <- args[5]
visualization_dir <- paste0(visualization_dir, output_stem)

############################
# Load in files
############################
sample_covariate_file <- paste0(processed_data_dir, "cell_covariates.txt")

############################
# Model Specification
############################
#model_stem <- paste0("eqtl_factorization_standard_eqtl_10.0_none_zscore_capped_eqtl_factorization_vi_ard_results_k_init_10_lambda_v_1_seed_2_var_param_1e-3_ratio_variance_std_True_permute_False_temper_")
eqtl_factorization_loading_file <- paste0(eqtl_results_dir, model_stem, "U_S.txt")


eqtl_factorization_factor_file <- paste0(eqtl_results_dir, model_stem, "V.txt")
pve_file <- paste0(eqtl_results_dir, model_stem, "factor_pve.txt")

pve <- as.numeric(read.table(pve_file, header=FALSE, sep="\t")$V1)

ordering <- order(pve, decreasing=TRUE)
#ordering <- ordering[1:3]
#print(ordering)



# Load in data
covariates <- read.table(sample_covariate_file, header=TRUE, sep="\t", comment.char="*")
loadings <- read.table(eqtl_factorization_loading_file, header=FALSE)
#factors <- read.table(eqtl_factorization_factor_file, header=FALSE)




loadings <- loadings[, ordering]
ordered_pve <- pve[ordering]

#factors <- factors[ordering,]

#######################################
# Covariate-loading heatmap
#######################################
output_file <- paste0(visualization_dir, "covariate_loading_correlation_heatmap.pdf")
heatmap <- make_covariate_loading_correlation_heatmap(covariates, loadings)
ggsave(heatmap, file=output_file, width=7.2, height=5.5, units="in")


#######################################
# PVE plot showing fraction of eqtl variance explained through each factor
#######################################
output_file <- paste0(visualization_dir, "fraction_of_eqtl_variance_explained_lineplot.pdf")
pve_plot <- make_pc_variance_explained_line_plot(ordered_pve)
ggsave(pve_plot, file=output_file, width=7.2, height=5.5, units="in")

#######################################
# Make histogram showing distribution of factor values for each factor
#######################################
output_file <- paste0(visualization_dir, "factor_distribution_histograms.pdf")
#hist <- make_factor_distribution_histograms(factors)
#ggsave(hist, file=output_file, width=7.2, height=7.5, units="in")

######################################
# Make loading boxplot with row for every factor colored by batch
#######################################
output_file <- paste0(visualization_dir, "loading_boxplot_with_row_for_every_factor_colored_by_plate.pdf")
boxplot <- make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(covariates$plate_id, loadings, "Plate ID")
ggsave(boxplot, file=output_file, width=10.2, height=20.5, units="in")


######################################
# Make loading vs loading scatter plot colored by real-valued covariate
#######################################
output_file <- paste0(visualization_dir, "loading1_vs_loading2_colored_by_g1_s_transition.pdf")
scatter <- make_loading_vs_loading_scatter_colored_by_real_valued_covariate(loadings[,1], loadings[,2], covariates$G1_S_transition, "Loading 1", "Loading 2", "G1-S transition")
ggsave(scatter, file=output_file, width=7.2, height=5.5, units="in")

######################################
# Make loading vs loading scatter plot colored by categorical covariate
#######################################
output_file <- paste0(visualization_dir, "loading1_vs_loading2_colored_by_day.pdf")
scatter <- make_loading_vs_loading_scatter_colored_by_categorical_covariate(loadings[,1], loadings[,2], covariates$day, "Loading 1", "Loading 2", "Differentiation day")
ggsave(scatter, file=output_file, width=7.2, height=5.5, units="in")

######################################
# Make loading boxplot with row for every factor colored by batch
#######################################
output_file <- paste0(visualization_dir, "loading_boxplot_with_row_for_every_factor_colored_by_day.pdf")
boxplot <- make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(covariates$day, loadings, "Plate ID")
ggsave(boxplot, file=output_file, width=10.2, height=20.5, units="in")

######################################
# Make loading boxplot colored by Differentiation day
#######################################
output_file <- paste0(visualization_dir, "loading_boxplot_colored_by_day.pdf")
boxplot <- make_loading_boxplot_plot_by_categorical_covariate(covariates$day, loadings, "Differentiation day")
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")


output_file <- paste0(visualization_dir, "loading_by_pseudotime_scatter_with_row_for_every_factor.pdf")
scatter_plot <-loading_expression_pc1_scatter_with_row_for_every_factor(covariates$PC1_top500hvgs, loadings[,1:2], "PC1")
ggsave(scatter_plot, file=output_file, width=7.2, height=10.5, units="in")




print('UMAP START')
umap_loadings = umap(loadings)$layout
saveRDS( umap_loadings, "umap_loadings.rds")
#umap_loadings <- readRDS("umap_loadings.rds")
print('UMAP DONE')


######################################
# Visualize UMAP scatter plot colored by individual
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_donor.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$donor, umap_loadings, "Known individual")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by cell type
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_day.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$day, umap_loadings, "day")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by cell type
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_experiment.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$experiment, umap_loadings, "experiment")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by cell type
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_plate.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$plate_id, umap_loadings, "plate")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


######################################
# Visualize UMAP scatter plot colored by number of cells
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_pc1_top500hvgs.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$PC1_top500hvgs, umap_loadings, "PC1")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")





