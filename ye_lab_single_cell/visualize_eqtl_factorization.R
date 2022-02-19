args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(cowplot)
library(umap)
library(ggplot2)
library(RColorBrewer)
library(sigmoid)
library(lme4)
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


make_umap_loading_scatter_plot_colored_by_categorical_variable <- function(covariates, umap_loadings, covariate_name) {
	df <- data.frame(loading_1=umap_loadings[,1], loading_2=umap_loadings[,2], covariate=factor(covariates))
	plotter <- ggplot(df) + 
	           geom_point(aes(x=loading_1, y=loading_2, color=covariate), size=.001) +
	           gtex_v8_figure_theme() + 
	           guides(colour = guide_legend(override.aes = list(size=2))) +
	           labs(x="SURGE UMAP 1", y = "SURGE UMAP 2", color=covariate_name) + 
	           guides(colour=guide_legend(nrow=3,byrow=TRUE, override.aes = list(size=2))) +
	           theme(legend.position="bottom") + 
	           theme(legend.text = element_text(size=8), legend.title = element_text(size=8))
	return(plotter)
}


make_umap_loading_scatter_plot_colored_by_real_valued_variable <- function(covariates, umap_loadings, covariate_name,point_size=.01) {
	df <- data.frame(loading_1=umap_loadings[,1], loading_2=umap_loadings[,2], covariate=covariates)
	plotter <- ggplot(df) + 
	           geom_point(aes(x=loading_1, y=loading_2, color=covariate), size=point_size) +
	           gtex_v8_figure_theme() + 
	           labs(x="SURGE UMAP 1", y = "SURGE UMAP 2", color=covariate_name) + 
             scale_colour_gradient2() +
	           scale_color_gradient(low="pink",high="blue") +
	           theme(legend.position="bottom") + 
	           theme(legend.text = element_text(size=8), legend.title = element_text(size=8))
	return(plotter)
}

######################################
# Make correlation heatmap correlating covariates with loadings
#######################################
make_pseudobulk_pc_loading_correlation_heatmap <- function(expr_pcs, loadings) {
  loadings <- as.matrix(loadings)
  expr_pcs <- as.matrix(expr_pcs)[,1:40]



  # Initialize PVE heatmap
  factor_colnames <- paste0("Factor", 1:(dim(loadings)[2]))
  factor_rownames <- paste0("Expression_pc", 1:(dim(expr_pcs)[2]))
  pve_map <- matrix(0, dim(expr_pcs)[2], dim(loadings)[2])
  colnames(pve_map) <- factor_colnames
  rownames(pve_map) <- factor_rownames


    # Loop through each PC, COV Pair and take correlation
    num_pcs <- dim(loadings)[2]
    num_expr_pcs <- dim(expr_pcs)[2]
    for (num_pc in 1:num_pcs) {
        for (num_expr_pc in 1:num_expr_pcs) {
            pc_vec <- loadings[,num_pc]
            expr_pc_vec <- expr_pcs[,num_expr_pc]
            #print(paste0(num_pc, " - ", num_cov))
            lin_model <- lm(pc_vec ~ expr_pc_vec)
            pve_map[num_expr_pc, num_pc] <- summary(lin_model)$adj.r.squared
        }
    }
    
    ord <- hclust( dist(scale(pve_map), method = "euclidean"), method = "ward.D" )$order

    melted_mat <- melt(pve_map)
    colnames(melted_mat) <- c("Expression_PC", "SURGE_PC","PVE")

    melted_mat$Expression_PC = factor(melted_mat$Expression_PC, levels=factor_rownames)
    melted_mat$SURGE_PC = factor(melted_mat$SURGE_PC, levels=factor_colnames)

    
    #levels(melted_mat$PC) = paste0("PC", 1:(length(levels(melted_mat$PC))))
    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=Expression_PC, y=SURGE_PC)) + geom_tile(aes(fill=PVE)) + scale_fill_gradient2(midpoint=-.05, guide="colorbar")
    heatmap <- heatmap + labs(y="",fill="VE")
    heatmap <- heatmap + theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
    # Save File
    return(heatmap)
}



######################################
# Make correlation heatmap correlating covariates with loadings
#######################################
make_pseudobulk_covariate_loading_correlation_heatmap <- function(covariates, loadings) {
  loadings <- as.matrix(loadings)


  valid_covariates <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33)
  covariate_type <- c("cat", "num", "cat", "cat", "cat", "cat", "num", "cat", "cat", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num", "num")

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
                ylim(0,max(variance_explained) + .002) + 
                scale_x_continuous(breaks=seq(0,(num_pcs),1)) +
                labs(x = "Latent context", y = "PVE") + 
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

generate_isg_signature_vector <- function(expr, gene_names) {
  marker_genes <- c("CD74", "IL2RB", "IRF8", "CD69", "MARCH1", "IFIT1")
  marker_gene_indices <- c()
  for (marker_gene_iter in 1:length(marker_genes)) {
    marker_gene = marker_genes[marker_gene_iter]
    if (marker_gene %in% marker_genes) {
       marker_gene_index <- which(gene_names==marker_gene)
       marker_gene_indices <- c(marker_gene_indices, marker_gene_index)
    }
  }

  base_expr <- expr[, marker_gene_indices[2]]

  for (iter in 1:length(marker_gene_indices)) {
    print(cor(base_expr, expr[, marker_gene_indices[iter]]))
  }


}

sldsc_enrichment_se_plot_over_continuous_domain <- function(df, trait_name, component_num, static_eqtl_trait_subset) {
  df$enrichment_lb <- df$enrichment - df$enrichment_std_err
  df$enrichment_ub <- df$enrichment + df$enrichment_std_err

  p <- ggplot(df,aes(x=component_position,y=enrichment)) + 
      geom_ribbon(aes(x=component_position,ymin=enrichment_lb,ymax=enrichment_ub),fill='thistle2')+
      geom_line(col='orchid3') + 
      gtex_v8_figure_theme() +
      labs(x=paste0("SURGE context ", component_num), y=paste0(trait_name,"\nS-LDSC enrichment")) + 
      geom_hline(yintercept=static_eqtl_trait_subset$enrichment[1], col='black', size=.5) + 
      geom_hline(yintercept=static_eqtl_trait_subset$enrichment[1] - static_eqtl_trait_subset$enrichment_std_err, col='black', linetype="dashed", size=.5) +
      geom_hline(yintercept=static_eqtl_trait_subset$enrichment[1] + static_eqtl_trait_subset$enrichment_std_err, col='black', linetype="dashed", size=.5)

  return(p)

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
per_cell_sldsc_results_dir <- args[6]
per_cell_3_component_sldsc_results_dir <- args[7]
component_gridspace_sldsc_results_dir <- args[8]
static_eqtl_sldsc_results_dir <- args[9]

############################
# Load in files
############################
sample_covariate_file <- paste0(processed_data_dir, "pseudobulk_scran_normalization_hvg_6000_regress_batch_True_individual_clustering_leiden_resolution_10.0_no_cap_15_sample_covariates.txt")
gene_names_file <- paste0(processed_data_dir, "pseudobulk_scran_normalization_hvg_6000_regress_batch_True_individual_clustering_leiden_resolution_10.0_no_cap_15_gene_names.txt")
gene_expr_file <- paste0(processed_data_dir, "pseudobulk_scran_normalization_hvg_6000_regress_batch_True_individual_clustering_leiden_resolution_10.0_no_cap_15_none_sample_norm_zscore_gene_norm_normalized_expression.txt")
gene_expr_pc_file <- paste0(processed_data_dir, "pseudobulk_scran_normalization_hvg_6000_regress_batch_True_individual_clustering_leiden_resolution_10.0_no_cap_15_none_sample_norm_zscore_gene_norm_pca_scores.txt")

if (FALSE) {
per_cell_sldsc_results_file <- paste0(per_cell_sldsc_results_dir, "per_cell_sldsc_results.txt")
per_cell_sldsc_blood_ma_results_file <- paste0(per_cell_sldsc_results_dir, "per_cell_sldsc_Blood_meta_analysis_results.txt")
per_cell_sldsc_immune_ma_results_file <- paste0(per_cell_sldsc_results_dir, "per_cell_sldsc_Immune_meta_analysis_results.txt")
per_cell_sldsc_non_blood_immune_ma_results_file <- paste0(per_cell_sldsc_results_dir, "per_cell_sldsc_Non_blood_immune_meta_analysis_results.txt")

per_cell_3_component_sldsc_results_file <- paste0(per_cell_3_component_sldsc_results_dir, "per_cell_sldsc_results.txt")
per_cell_3_component_sldsc_blood_ma_results_file <- paste0(per_cell_3_component_sldsc_results_dir, "per_cell_sldsc_Blood_meta_analysis_results.txt")
per_cell_3_component_sldsc_immune_ma_results_file <- paste0(per_cell_3_component_sldsc_results_dir, "per_cell_sldsc_Immune_meta_analysis_results.txt")
per_cell_3_component_sldsc_non_blood_immune_ma_results_file <- paste0(per_cell_3_component_sldsc_results_dir, "per_cell_sldsc_Non_blood_immune_meta_analysis_results.txt")

component_gridspace_sldsc_results_file <- paste0(component_gridspace_sldsc_results_dir, "component_gridspace_sldsc_results.txt")

static_eqtl_sldsc_results_file <- paste0(static_eqtl_sldsc_results_dir, "static_eqtl_sldsc_results.txt")
}


############################
# Model Specification
############################
#model_stem <- paste0("eqtl_factorization_standard_eqtl_10.0_none_zscore_capped_eqtl_factorization_vi_ard_results_k_init_10_lambda_v_1_seed_2_var_param_1e-3_ratio_variance_std_True_permute_False_temper_")
eqtl_factorization_loading_file <- paste0(eqtl_results_dir, model_stem, "theta_normalized.txt")
eqtl_factorization_loading_file <- paste0(eqtl_results_dir, model_stem, "U_S.txt")


eqtl_factorization_factor_file <- paste0(eqtl_results_dir, model_stem, "V.txt")
pve_file <- paste0(eqtl_results_dir, model_stem, "factor_pve.txt")

pve <- as.numeric(read.table(pve_file, header=FALSE, sep="\t")$V1)

ordering <- order(pve, decreasing=TRUE)
#ordering <- ordering[1:4]
#print(ordering)



# Load in data
covariates <- read.table(sample_covariate_file, header=TRUE, sep="\t")

if (FALSE) {
per_cell_sldsc_results <- read.table(per_cell_sldsc_results_file, header=TRUE, sep="\t")
per_cell_blood_ma_sldsc_results <- read.table(per_cell_sldsc_blood_ma_results_file, header=TRUE, sep="\t")
per_cell_immune_ma_sldsc_results <- read.table(per_cell_sldsc_immune_ma_results_file, header=TRUE, sep="\t")
per_cell_non_blood_immune_ma_sldsc_results <- read.table(per_cell_sldsc_non_blood_immune_ma_results_file, header=TRUE, sep="\t")

per_cell_3_component_sldsc_results <- read.table(per_cell_3_component_sldsc_results_file, header=TRUE, sep="\t")
per_cell_3_component_blood_ma_sldsc_results <- read.table(per_cell_3_component_sldsc_blood_ma_results_file, header=TRUE, sep="\t")
per_cell_3_component_immune_ma_sldsc_results <- read.table(per_cell_3_component_sldsc_immune_ma_results_file, header=TRUE, sep="\t")
per_cell_3_component_non_blood_immune_ma_sldsc_results <- read.table(per_cell_3_component_sldsc_non_blood_immune_ma_results_file, header=TRUE, sep="\t")


component_gridspace_sldsc_results <- read.table(component_gridspace_sldsc_results_file, header=TRUE, sep="\t")
static_eqtl_sldsc_results <- read.table(static_eqtl_sldsc_results_file, header=TRUE, sep="\t")
}

# Change "nan" to "monocyte"
covariates$ct_cov_mode = as.character(covariates$ct_cov_mode)
covariates$ct_cov_mode[covariates$ct_cov_mode == "nan"] = rep("monocyte", sum(covariates$ct_cov_mode == "nan"))
covariates$ct_cov_mode = factor(covariates$ct_cov_mode)


loadings <- read.table(eqtl_factorization_loading_file, header=FALSE)
#factors <- read.table(eqtl_factorization_factor_file, header=FALSE)

gene_names <- read.table(gene_names_file, header=FALSE)$V1


expr <- read.table(gene_expr_file, header=FALSE)
saveRDS( expr, "expr.rds")
#expr <- readRDS("expr.rds")

expr_pcs <- read.table(gene_expr_pc_file, header=FALSE)
saveRDS(expr_pcs, "expr_pcs.rds")
#expr_pcs <- readRDS("expr_pcs.rds")


umap_expr_file <- paste0(processed_data_dir, "temp_15_umap.txt")
umap_expr <- read.table(umap_expr_file, header=FALSE, sep="\t")
#print("DONE")

#gene_expr_pc_file <- "/work-zfs/abattle4/bstrober/qtl_factorization/ye_lab_single_cell/eqtl_factorization_results/eqtl_factorization_standard_eqtl_hvg_6000_10.0_no_cap_15_none_zscore_factorization_vi_ard_results_k_init_30_seed_1_warmup_3000_ratio_variance_std_True_permute_False_lambda_1_round_geno_True_temper_U_S.txt"
#expr_pcs <- read.table(gene_expr_pc_file, header=FALSE)
#umap_expr = umap(expr_pcs)$layout
#saveRDS(umap_expr, "umap_expr_loadings.rds")


loadings <- loadings[, ordering]
ordered_pve <- pve[ordering]



covariates$ct_by_status = factor(paste0(covariates$ct_cov_mode, "_", covariates$SLE_status))
covariates$cg_by_status = factor(paste0(covariates$cg_cov_mode, "_", covariates$SLE_status))
covariates$ct_by_pop = factor(paste0(covariates$ct_cov_mode, "_", covariates$pop_cov))
covariates$cg_by_pop = factor(paste0(covariates$cg_cov_mode, "_", covariates$pop_cov))

#fit <- lmer(expr_pcs[,3] ~ loadings[,1] + loadings[,2] + loadings[,3] + loadings[,4] + loadings[,5] + loadings[,6] + loadings[,7] + loadings[,8] + loadings[,9] + loadings[,10] + (1|Z))
#fit2 <- lm(expr_pcs[,3] ~ predict(fit))
#print(summary(fit2))

#######################################
# Generate isg signature vector
#######################################
#isg_signature_vector = generate_isg_signature_vector(expr, gene_names)

#######################################
# Expression UMAP scatter colored by surge factor loadings
#######################################
surge_num <- 1
surge_vec <- loadings[,surge_num]
output_file <- paste0(visualization_dir, "expr_pc_umap_loading_scatter_colored_by_surge_factor_num_", surge_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(surge_vec, umap_expr, paste0("SURGE ", surge_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")
surge_num <- 2
surge_vec <- loadings[,surge_num]
output_file <- paste0(visualization_dir, "expr_pc_umap_loading_scatter_colored_by_surge_factor_num_", surge_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(surge_vec, umap_expr, paste0("SURGE ", surge_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")
surge_num <- 3
surge_vec <- loadings[,surge_num]
output_file <- paste0(visualization_dir, "expr_pc_umap_loading_scatter_colored_by_surge_factor_num_", surge_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(surge_vec, umap_expr, paste0("SURGE ", surge_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")
surge_num <- 4
surge_vec <- loadings[,surge_num]
output_file <- paste0(visualization_dir, "expr_pc_umap_loading_scatter_colored_by_surge_factor_num_", surge_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(surge_vec, umap_expr, paste0("SURGE ", surge_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")
surge_num <- 5
surge_vec <- loadings[,surge_num]
output_file <- paste0(visualization_dir, "expr_pc_umap_loading_scatter_colored_by_surge_factor_num_", surge_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(surge_vec, umap_expr, paste0("SURGE ", surge_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")
surge_num <- 6
surge_vec <- loadings[,surge_num]
output_file <- paste0(visualization_dir, "expr_pc_umap_loading_scatter_colored_by_surge_factor_num_", surge_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(surge_vec, umap_expr, paste0("SURGE ", surge_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")
surge_num <- 7
surge_vec <- loadings[,surge_num]
output_file <- paste0(visualization_dir, "expr_pc_umap_loading_scatter_colored_by_surge_factor_num_", surge_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(surge_vec, umap_expr, paste0("SURGE ", surge_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")
surge_num <- 8
surge_vec <- loadings[,surge_num]
output_file <- paste0(visualization_dir, "expr_pc_umap_loading_scatter_colored_by_surge_factor_num_", surge_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(surge_vec, umap_expr, paste0("SURGE ", surge_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")
surge_num <- 9
surge_vec <- loadings[,surge_num]
output_file <- paste0(visualization_dir, "expr_pc_umap_loading_scatter_colored_by_surge_factor_num_", surge_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(surge_vec, umap_expr, paste0("SURGE ", surge_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")
surge_num <- 10
surge_vec <- loadings[,surge_num]
output_file <- paste0(visualization_dir, "expr_pc_umap_loading_scatter_colored_by_surge_factor_num_", surge_num, ".pdf")#
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(surge_vec, umap_expr, paste0("SURGE ", surge_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


######################################
# Visualize UMAP scatter plot colored by cell type
#######################################
output_file <- paste0(visualization_dir, "expr_pc_umap_loading_scatter_colored_by_cell_type.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$cg_cov_mode, umap_expr, "cell type")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by cell type
#######################################
output_file <- paste0(visualization_dir, "expr_pc_umap_loading_scatter_colored_by_cell_type2.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$ct_cov_mode, umap_expr, "cell type")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by cell type
#######################################
output_file <- paste0(visualization_dir, "expr_pc_umap_loading_scatter_colored_by_individual.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$ind_cov, umap_expr, "cell type")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")



#######################################
# PVE plot showing fraction of eqtl variance explained through each factor
#######################################
output_file <- paste0(visualization_dir, "fraction_of_eqtl_variance_explained_lineplot.pdf")
pve_plot <- make_pc_variance_explained_line_plot(ordered_pve[1:3])
ggsave(pve_plot, file=output_file, width=7.2, height=5.5, units="in")

#######################################
# Make histogram showing distribution of factor values for each factor
#######################################
output_file <- paste0(visualization_dir, "factor_distribution_histograms.pdf")
#hist <- make_factor_distribution_histograms(factors)
#ggsave(hist, file=output_file, width=7.2, height=7.5, units="in")


loading_num <- 1
output_file <- paste0(visualization_dir, "histogram_of_loadings_", loading_num, "_for_each_cell_type_stratefied_by_sle_status.pdf")
histo <- make_histogram_of_loadings_for_each_cell_type_stratefied_by_sle_status(loadings[,loading_num], covariates$cg_cov_mode, covariates$SLE_status, loading_num)
ggsave(histo, file=output_file, width=7.2, height=12, units="in")

loading_num <- 1
cell_type <- "monocyte"
cell_type_indices = as.character(covariates$ct_cov_mode) == cell_type
loading_vec <- loadings[,loading_num]
output_file <- paste0(visualization_dir, "histogram_of_loadings_", loading_num, "_for_", cell_type, "_stratefied_by_sle_status.pdf")
histy <- make_histogram_of_loadings_for_cell_type_stratefied_by_sle_status(loading_vec[cell_type_indices], covariates$SLE_status[cell_type_indices], cell_type, loading_num)
ggsave(histy, file=output_file, width=7.2, height=4, units="in")


######################################
# Make loading boxplot colored by Ancestry
#######################################
output_file <- paste0(visualization_dir, "loading_boxplot_colored_by_ancestry.pdf")
boxplot <- make_loading_boxplot_plot_by_categorical_covariate(covariates$pop_cov, loadings, "Known ancestry")
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")

######################################
# Make loading boxplot colored by Ancestry
#######################################
output_file <- paste0(visualization_dir, "loading_boxplot_colored_by_sle_status.pdf")
boxplot <- make_loading_boxplot_plot_by_categorical_covariate(covariates$SLE_status, loadings, "SLE status")
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")


######################################
# Make loading boxplot colored by Ancestry
#######################################
output_file <- paste0(visualization_dir, "loading_boxplot_colored_by_cell_type.pdf")
boxplot <- make_loading_boxplot_plot_by_categorical_covariate(covariates$cg_cov_mode, loadings, "Cell Type")
ggsave(boxplot, file=output_file, width=12.2, height=5.5, units="in")

######################################
# Make loading boxplot colored by Ancestry
#######################################
output_file <- paste0(visualization_dir, "loading_boxplot_colored_by_cell_typeXstatus.pdf")
boxplot <- make_loading_boxplot_plot_by_categorical_covariate(covariates$ct_by_status, loadings, "Cell Type X status")
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")

######################################
# Make loading boxplot colored by Ancestry
#######################################
output_file <- paste0(visualization_dir, "loading_boxplot_colored_by_cell_type2Xstatus.pdf")
boxplot <- make_loading_boxplot_plot_by_categorical_covariate(covariates$cg_by_status, loadings, "Cell Type (fine grained) X status")
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")

######################################
# Make loading boxplot colored by Ancestry
#######################################
output_file <- paste0(visualization_dir, "loading_boxplot_colored_by_cell_typeXpop.pdf")
boxplot <- make_loading_boxplot_plot_by_categorical_covariate(covariates$ct_by_pop, loadings, "Cell Type X pop")
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")

######################################
# Make loading boxplot colored by Ancestry
#######################################
output_file <- paste0(visualization_dir, "loading_boxplot_colored_by_cell_type2Xpop.pdf")
boxplot <- make_loading_boxplot_plot_by_categorical_covariate(covariates$cg_by_pop, loadings, "Cell Type (fine grained) X pop")
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")



######################################
# Make loading boxplot colored by Ancestry
#######################################
output_file <- paste0(visualization_dir, "loading_boxplot_colored_by_cell_type_fine_res.pdf")
boxplot <- make_loading_boxplot_plot_by_categorical_covariate(covariates$ct_cov_mode, loadings, "Cell Type")
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
output_file <- paste0(visualization_dir, "eqtl_factorization_covariate_pve_heatmap.pdf")
ggsave(heatmap, file=output_file, width=7.2, height=10, units="in")


##########################
# PCA-covariate heatmap for pseudobulk data
##########################
heatmap <- make_pseudobulk_pc_loading_correlation_heatmap(expr_pcs, loadings)
output_file <- paste0(visualization_dir, "eqtl_factorization_pca_pve_heatmap.pdf")
ggsave(heatmap, file=output_file, width=7.2, height=10, units="in")


######################################
# Make loading boxplot with row for every factor colored by individual
#######################################
output_file <- paste0(visualization_dir, "loading_boxplot_with_row_for_every_factor_colored_by_individual.pdf")
boxplot <- make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(covariates$ind_cov, loadings, "Individual")
ggsave(boxplot, file=output_file, width=10.2, height=20.5, units="in")

######################################
# Make loading boxplot with row for every factor colored by batch
#######################################
output_file <- paste0(visualization_dir, "loading_boxplot_with_row_for_every_factor_colored_by_batch.pdf")
boxplot <- make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(covariates$batch_cov, loadings, "Batch")
ggsave(boxplot, file=output_file, width=10.2, height=20.5, units="in")


######################################
# Make loading boxplot with row for every factor colored by individual limited to B cells
#######################################
cell_indices <- as.character(covariates$cg_cov_mode) == "B"
output_file <- paste0(visualization_dir, "loading_boxplot_with_row_for_every_factor_colored_by_individual_for_only_B_cells.pdf")
boxplot <- make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(covariates$ind_cov[cell_indices], loadings[cell_indices,], "Individual")
ggsave(boxplot, file=output_file, width=10.2, height=20.5, units="in")

######################################
# Make loading boxplot with row for every factor colored by batch limited to B cells
#######################################
output_file <- paste0(visualization_dir, "loading_boxplot_with_row_for_every_factor_colored_by_batch_for_only_B_cells.pdf")
boxplot <- make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(covariates$batch_cov[cell_indices], loadings[cell_indices,], "Batch")
ggsave(boxplot, file=output_file, width=10.2, height=20.5, units="in")


######################################
# Make loading boxplot with row for every factor colored by individual limited to NK cells
#######################################
cell_indices <- as.character(covariates$cg_cov_mode) == "NK"
output_file <- paste0(visualization_dir, "loading_boxplot_with_row_for_every_factor_colored_by_individual_for_only_NK_cells.pdf")
boxplot <- make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(covariates$ind_cov[cell_indices], loadings[cell_indices,], "Individual")
ggsave(boxplot, file=output_file, width=10.2, height=20.5, units="in")

######################################
# Make loading boxplot with row for every factor colored by batch limited to NK cells
#######################################
output_file <- paste0(visualization_dir, "loading_boxplot_with_row_for_every_factor_colored_by_batch_for_only_NK_cells.pdf")
boxplot <- make_loading_boxplot_plot_with_row_for_every_factor_by_categorical_covariate(covariates$batch_cov[cell_indices], loadings[cell_indices,], "Batch")
ggsave(boxplot, file=output_file, width=10.2, height=20.5, units="in")




if (FALSE) {
trait_arr <- unique(as.character(component_gridspace_sldsc_results$trait_name))


for (trait_iter in 1:length(trait_arr)) {
  trait_name <- trait_arr[trait_iter]

  output_file <- paste0(visualization_dir, "continous_", trait_name, "_sldsc_enrichment_plot.pdf")

  static_eqtl_trait_subset <- static_eqtl_sldsc_results[as.character(static_eqtl_sldsc_results$trait_name) == trait_name,]

  component_num <- 1
  indices <- (component_gridspace_sldsc_results$trait_name == trait_name) & (component_gridspace_sldsc_results$component_num == component_num)
  continous_sldsc_enrichmennt_plot1 <- sldsc_enrichment_se_plot_over_continuous_domain(component_gridspace_sldsc_results[indices,], trait_name, component_num, static_eqtl_trait_subset)

  component_num <- 2
  indices <- (component_gridspace_sldsc_results$trait_name == trait_name) & (component_gridspace_sldsc_results$component_num == component_num)
  continous_sldsc_enrichmennt_plot2 <- sldsc_enrichment_se_plot_over_continuous_domain(component_gridspace_sldsc_results[indices,], trait_name, component_num, static_eqtl_trait_subset)

  component_num <- 3
  indices <- (component_gridspace_sldsc_results$trait_name == trait_name) & (component_gridspace_sldsc_results$component_num == component_num)
  continous_sldsc_enrichmennt_plot3 <- sldsc_enrichment_se_plot_over_continuous_domain(component_gridspace_sldsc_results[indices,], trait_name, component_num, static_eqtl_trait_subset)

  continous_sldsc_enrichmennt_plot_joint <- plot_grid(continous_sldsc_enrichmennt_plot1,continous_sldsc_enrichmennt_plot2, continous_sldsc_enrichmennt_plot3, ncol=1)

  ggsave(continous_sldsc_enrichmennt_plot_joint, file=output_file, width=7.2, height=8.0, units="in")

}



output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_sldsc_blood_meta_enrichment_3_component.pdf")
trait_enrichment <- per_cell_3_component_blood_ma_sldsc_results$enrichment
indices = !is.na(trait_enrichment)
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(trait_enrichment[indices], umap_loadings_3_comp[indices,], paste0("S-LDSC enrichment (blood meta analysis)"), point_size=1)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_sldsc_blood_meta_tau_3_component.pdf")
trait_enrichment <- per_cell_3_component_blood_ma_sldsc_results$tau_star
indices = !is.na(trait_enrichment)
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(trait_enrichment[indices], umap_loadings_3_comp[indices,], paste0("S-LDSC Tau (blood meta analysis)"), point_size=1)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_sldsc_immune_meta_enrichment_3_component.pdf")
trait_enrichment <- per_cell_3_component_immune_ma_sldsc_results$enrichment
indices = !is.na(trait_enrichment)
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(trait_enrichment[indices], umap_loadings_3_comp[indices,], paste0("S-LDSC enrichment (immune meta analysis)"), point_size=1)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_sldsc_immune_meta_tau_3_component.pdf")
trait_enrichment <- per_cell_3_component_immune_ma_sldsc_results$tau_star
indices = !is.na(trait_enrichment)
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(trait_enrichment[indices], umap_loadings_3_comp[indices,], paste0("S-LDSC Tau (immune meta analysis)"), point_size=1)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_sldsc_non_blood_immune_meta_enrichment_3_component.pdf")
trait_enrichment <- per_cell_3_component_non_blood_immune_ma_sldsc_results$enrichment
indices = !is.na(trait_enrichment)
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(trait_enrichment[indices], umap_loadings_3_comp[indices,], paste0("S-LDSC enrichment (non-blood-immune meta analysis)"), point_size=1)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_sldsc_non_blood_immune_meta_tau_3_component.pdf")
trait_enrichment <- per_cell_3_component_non_blood_immune_ma_sldsc_results$tau_star
indices = !is.na(trait_enrichment)
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(trait_enrichment[indices], umap_loadings_3_comp[indices,], paste0("S-LDSC Tau (non-blood-immune immune meta analysis)"), point_size=1)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")





trait_names = c('ukbb_blood_monocyte_count', 'ukbb_blood_lymphocyte_count', 'ukbb_bmi', 'ukbb_eczema', 'ukbb_blood_eosinophil_count', 'ukbb_blood_high_light_scatter_reticulotye_count', 'ukbb_blood_mean_corpuscular_hemoglobin', 'ukbb_blood_platelet_vol', 'ukbb_blood_platelet_count', 'ukbb_blood_red_count', 'ukbb_blood_white_count', 'ukbb_height', 'ukbb_T2D', 'Celiac', 'Crohns', 'Ulcerative_Colitis', 'Rheumatoid_Arthritis', 'Lupus', 'IBD', 'Multiple_sclerosis', 'PBC', 'CAD', 'Bipolar', 'Alzheimer', 'Schizophrenia')

for (trait_index in 1:length(trait_names)) {
  trait_name <- trait_names[trait_index]

  output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_sldsc_", trait_name, "_enrichment_3_component.pdf")
  trait_enrichment <- per_cell_3_component_sldsc_results[paste0(trait_name, "_enrichment")][,1]
  indices = !is.na(trait_enrichment)
  umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(trait_enrichment[indices], umap_loadings_3_comp[indices,], paste0("S-LDSC enrichment (", trait_name, ")"), point_size=1)
  ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

  output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_sldsc_", trait_name, "_tau_3_component.pdf")
  trait_tau <- per_cell_3_component_sldsc_results[paste0(trait_name, "_tau_star")][,1]
  indices = !is.na(trait_tau)
  umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(trait_tau[indices], umap_loadings_3_comp[indices,], paste0("S-LDSC Tau (", trait_name, ")"), point_size=1)
  ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


}

}

print('UMAP START')
umap_loadings = umap(loadings)$layout
saveRDS( umap_loadings, "umap_loadings.rds")
#umap_loadings <- readRDS("umap_loadings.rds")
print('UMAP DONE')

#indices <- (abs(umap_loadings[,1]) < 8) & (abs(umap_loadings[,2]) < 8)
#umap_loadings <- umap_loadings[indices,]
#covariates <- covariates[indices,]
#loadings <- loadings[indices,]
#expr_pcs <- expr_pcs[indices,]
#expr <- expr[indices,]

if (FALSE) {
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_sldsc_blood_meta_enrichment.pdf")
trait_enrichment <- per_cell_blood_ma_sldsc_results$enrichment
indices = !is.na(trait_enrichment)
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(trait_enrichment[indices], umap_loadings[indices,], paste0("S-LDSC enrichment (blood meta analysis)"), point_size=1)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_sldsc_blood_meta_tau.pdf")
trait_enrichment <- per_cell_blood_ma_sldsc_results$tau_star
indices = !is.na(trait_enrichment)
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(trait_enrichment[indices], umap_loadings[indices,], paste0("S-LDSC Tau (blood meta analysis)"), point_size=1)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_sldsc_immune_meta_enrichment.pdf")
trait_enrichment <- per_cell_immune_ma_sldsc_results$enrichment
indices = !is.na(trait_enrichment)
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(trait_enrichment[indices], umap_loadings[indices,], paste0("S-LDSC enrichment (immune meta analysis)"), point_size=1)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_sldsc_immune_meta_tau.pdf")
trait_enrichment <- per_cell_immune_ma_sldsc_results$tau_star
indices = !is.na(trait_enrichment)
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(trait_enrichment[indices], umap_loadings[indices,], paste0("S-LDSC Tau (immune meta analysis)"), point_size=1)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_sldsc_non_blood_immune_meta_enrichment.pdf")
trait_enrichment <- per_cell_non_blood_immune_ma_sldsc_results$enrichment
indices = !is.na(trait_enrichment)
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(trait_enrichment[indices], umap_loadings[indices,], paste0("S-LDSC enrichment (non-blood-immune meta analysis)"), point_size=1)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_sldsc_non_blood_immune_meta_tau.pdf")
trait_enrichment <- per_cell_non_blood_immune_ma_sldsc_results$tau_star
indices = !is.na(trait_enrichment)
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(trait_enrichment[indices], umap_loadings[indices,], paste0("S-LDSC Tau (non-blood-immune immune meta analysis)"), point_size=1)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")





trait_names = c('ukbb_blood_monocyte_count', 'ukbb_blood_lymphocyte_count', 'ukbb_bmi', 'ukbb_eczema', 'ukbb_blood_eosinophil_count', 'ukbb_blood_high_light_scatter_reticulotye_count', 'ukbb_blood_mean_corpuscular_hemoglobin', 'ukbb_blood_platelet_vol', 'ukbb_blood_platelet_count', 'ukbb_blood_red_count', 'ukbb_blood_white_count', 'ukbb_height', 'ukbb_T2D', 'Celiac', 'Crohns', 'Ulcerative_Colitis', 'Rheumatoid_Arthritis', 'Lupus', 'IBD', 'Multiple_sclerosis', 'PBC', 'CAD', 'Bipolar', 'Alzheimer', 'Schizophrenia')

for (trait_index in 1:length(trait_names)) {
  trait_name <- trait_names[trait_index]

  output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_sldsc_", trait_name, "_enrichment.pdf")
  trait_enrichment <- per_cell_sldsc_results[paste0(trait_name, "_enrichment")][,1]
  indices = !is.na(trait_enrichment)
  umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(trait_enrichment[indices], umap_loadings[indices,], paste0("S-LDSC enrichment (", trait_name, ")"), point_size=1)
  ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

  output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_sldsc_", trait_name, "_tau.pdf")
  trait_tau <- per_cell_sldsc_results[paste0(trait_name, "_tau_star")][,1]
  indices = !is.na(trait_tau)
  umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(trait_tau[indices], umap_loadings[indices,], paste0("S-LDSC Tau (", trait_name, ")"), point_size=1)
  ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


}
}


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
# Visualize UMAP scatter plot colored by cell type
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_cell_type2.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$ct_cov_mode, umap_loadings, "cell type")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")



######################################
# Visualize UMAP scatter plot colored by cell type
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_status.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$Status, umap_loadings, "Status")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by cell type
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_sle_status.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$SLE_status, umap_loadings, "SLE Status")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


######################################
# Visualize UMAP scatter plot colored by number of cells
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_log_donor_isg_score.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(log(covariates$donor_isg_score), umap_loadings, "log(Donor ISG score)")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by number of cells
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_donor_isg_score.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable((covariates$donor_isg_score), umap_loadings, "Donor ISG score")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by number of cells
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_num_cells.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$num_cells, umap_loadings, "Number of cells")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by number of cells
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_T8_fraction.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$T8_fraction, umap_loadings, "T8 fraction")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by number of cells
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_T4_fraction.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$T4_fraction, umap_loadings, "T4 fraction")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by number of cells
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_NK_fraction.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$NK_fraction, umap_loadings, "NK fraction")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by number of cells
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_cM_fraction.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$cM_fraction, umap_loadings, "cM fraction")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by number of cells
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_B_fraction.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(covariates$B_fraction, umap_loadings, "B fraction")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by number of cells
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_pop_cov.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$pop_cov, umap_loadings, "pop cov")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

######################################
# Visualize UMAP scatter plot colored by number of cells
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_sex.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$Sex, umap_loadings, "sex")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


######################################
# Visualize UMAP scatter plot colored by number of cells
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_batch.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$batch_cov, umap_loadings, "Batch")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


######################################
# Visualize UMAP scatter plot colored by number of cells
#######################################
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_ancestry.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_categorical_variable(covariates$pop_cov, umap_loadings, "Known Ancestry")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

print('here')
######################################
# Visualize UMAP scatter plot based on latent factors
#######################################
lf_num <- 1
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_ef_lf_", lf_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(sigmoid(loadings[,lf_num])-.5, umap_loadings, paste0("SURGE Latent context ", lf_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

lf_num <- 2
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_ef_lf_", lf_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(sigmoid(loadings[,lf_num])-.5, umap_loadings, paste0("SURGE Latent context ", lf_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

lf_num <- 3
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_ef_lf_", lf_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(sigmoid(loadings[,lf_num])-.5, umap_loadings, paste0("SURGE Latent context ", lf_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

lf_num <- 4
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_ef_lf_", lf_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(sigmoid(loadings[,lf_num])-.5, umap_loadings, paste0("SURGE Latent context ", lf_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

lf_num <- 5
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_ef_lf_", lf_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(sigmoid(loadings[,lf_num])-.5, umap_loadings, paste0("SURGE Latent context ", lf_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

lf_num <- 6
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_ef_lf_", lf_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(sigmoid(loadings[,lf_num])-.5, umap_loadings, paste0("SURGE Latent context ", lf_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


lf_num <- 7
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_ef_lf_", lf_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(sigmoid(loadings[,lf_num])-.5, umap_loadings, paste0("SURGE Latent context ", lf_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


lf_num <- 8
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_ef_lf_", lf_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(sigmoid(loadings[,lf_num])-.5, umap_loadings, paste0("SURGE Latent context ", lf_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


lf_num <- 9
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_ef_lf_", lf_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(sigmoid(loadings[,lf_num])-.5, umap_loadings, paste0("SURGE Latent context ", lf_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

lf_num <- 10
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_ef_lf_", lf_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(sigmoid(loadings[,lf_num])-.5, umap_loadings, paste0("SURGE Latent context ", lf_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")



######################################
# Visualize UMAP scatter plot based on expression pcs
#######################################
expression_pc_num <- 1
expr_pc_vec <- expr_pcs[,expression_pc_num]
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_expression_pc_num_", expression_pc_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(expr_pc_vec, umap_loadings, paste0("Expression PC ", expression_pc_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

expression_pc_num <- 2
expr_pc_vec <- expr_pcs[,expression_pc_num]
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_expression_pc_num_", expression_pc_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(expr_pc_vec, umap_loadings, paste0("Expression PC ", expression_pc_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

expression_pc_num <- 3
expr_pc_vec <- expr_pcs[,expression_pc_num]
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_expression_pc_num_", expression_pc_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(expr_pc_vec, umap_loadings, paste0("Expression PC ", expression_pc_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

expression_pc_num <- 4
expr_pc_vec <- expr_pcs[,expression_pc_num]
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_expression_pc_num_", expression_pc_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(expr_pc_vec, umap_loadings, paste0("Expression PC ", expression_pc_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

expression_pc_num <- 5
expr_pc_vec <- expr_pcs[,expression_pc_num]
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_expression_pc_num_", expression_pc_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(expr_pc_vec, umap_loadings, paste0("Expression PC ", expression_pc_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

expression_pc_num <- 6
expr_pc_vec <- expr_pcs[,expression_pc_num]
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_expression_pc_num_", expression_pc_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(expr_pc_vec, umap_loadings, paste0("Expression PC ", expression_pc_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

expression_pc_num <- 7
expr_pc_vec <- expr_pcs[,expression_pc_num]
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_expression_pc_num_", expression_pc_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(expr_pc_vec, umap_loadings, paste0("Expression PC ", expression_pc_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

expression_pc_num <- 8
expr_pc_vec <- expr_pcs[,expression_pc_num]
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_expression_pc_num_", expression_pc_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(expr_pc_vec, umap_loadings, paste0("Expression PC ", expression_pc_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


expression_pc_num <- 9
expr_pc_vec <- expr_pcs[,expression_pc_num]
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_expression_pc_num_", expression_pc_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(expr_pc_vec, umap_loadings, paste0("Expression PC ", expression_pc_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

expression_pc_num <- 10
expr_pc_vec <- expr_pcs[,expression_pc_num]
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_expression_pc_num_", expression_pc_num, ".pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(expr_pc_vec, umap_loadings, paste0("Expression PC ", expression_pc_num))
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_avg_cell_isg_score.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable((covariates$avg_cell_isg_score), umap_loadings, "avg cell ISG score")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_log_avg_cell_isg_score.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(log(covariates$avg_cell_isg_score - min(covariates$avg_cell_isg_score) + 1), umap_loadings, "log(avg cell ISG score)")
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")
######################################
# Visualize UMAP scatter plot based on marker genes
#######################################


marker_gene = "ISG15"
marker_gene_index <- which(gene_names==marker_gene)
marker_gene_expr_vec <- expr[, marker_gene_index]
marker_gene_expr_vec[marker_gene_expr_vec > 4.0] = 4.0
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_", marker_gene, "_expression.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(marker_gene_expr_vec, umap_loadings, marker_gene)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


marker_gene = "CD40LG"
marker_gene_index <- which(gene_names==marker_gene)
marker_gene_expr_vec <- expr[, marker_gene_index]
marker_gene_expr_vec[marker_gene_expr_vec > 4.0] = 4.0
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_", marker_gene, "_expression.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(marker_gene_expr_vec, umap_loadings, marker_gene)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

marker_gene = "CD8A"
marker_gene_index <- which(gene_names==marker_gene)
marker_gene_expr_vec <- expr[, marker_gene_index]
marker_gene_expr_vec[marker_gene_expr_vec > 4.0] = 4.0
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_", marker_gene, "_expression.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(marker_gene_expr_vec, umap_loadings, marker_gene)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

marker_gene = "PRF1"
marker_gene_index <- which(gene_names==marker_gene)
marker_gene_expr_vec <- expr[, marker_gene_index]
marker_gene_expr_vec[marker_gene_expr_vec > 6.0] = 6.0
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_", marker_gene, "_expression.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(marker_gene_expr_vec, umap_loadings, marker_gene)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


marker_gene = "BANK1"
marker_gene_index <- which(gene_names==marker_gene)
marker_gene_expr_vec <- expr[, marker_gene_index]
marker_gene_expr_vec[marker_gene_expr_vec > 6.0] = 6.0
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_", marker_gene, "_expression.pdf")
umap_scatter_bank1 <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(marker_gene_expr_vec, umap_loadings, marker_gene)
ggsave(umap_scatter_bank1, file=output_file, width=7.2, height=6.0, units="in")


marker_gene = "RTKN2"
marker_gene_index <- which(gene_names==marker_gene)
marker_gene_expr_vec <- expr[, marker_gene_index]
marker_gene_expr_vec[marker_gene_expr_vec > 4.0] = 4.0
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_", marker_gene, "_expression.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(marker_gene_expr_vec, umap_loadings, marker_gene)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

marker_gene = "CD8B"
marker_gene_index <- which(gene_names==marker_gene)
marker_gene_expr_vec <- expr[, marker_gene_index]
marker_gene_expr_vec[marker_gene_expr_vec > 4.0] = 4.0
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_", marker_gene, "_expression.pdf")
umap_scatter_cd8b <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(marker_gene_expr_vec, umap_loadings, marker_gene)
ggsave(umap_scatter_cd8b, file=output_file, width=7.2, height=6.0, units="in")

marker_gene = "ISG15"
marker_gene_index <- which(gene_names==marker_gene)
marker_gene_expr_vec <- expr[, marker_gene_index]
marker_gene_expr_vec[marker_gene_expr_vec > 4.0] = 4.0
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_", marker_gene, "_expression.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(marker_gene_expr_vec, umap_loadings, marker_gene)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

marker_gene = "IFNG"
marker_gene_index <- which(gene_names==marker_gene)
marker_gene_expr_vec <- expr[, marker_gene_index]
marker_gene_expr_vec[marker_gene_expr_vec > 4.0] = 4.0
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_", marker_gene, "_expression.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(marker_gene_expr_vec, umap_loadings, marker_gene)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")


marker_gene = "MX1"
marker_gene_index <- which(gene_names==marker_gene)
marker_gene_expr_vec <- expr[, marker_gene_index]
marker_gene_expr_vec[marker_gene_expr_vec > 4.0] = 4.0
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_", marker_gene, "_expression.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(marker_gene_expr_vec, umap_loadings, marker_gene)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

marker_gene = "LAG3"
marker_gene_index <- which(gene_names==marker_gene)
marker_gene_expr_vec <- expr[, marker_gene_index]
marker_gene_expr_vec[marker_gene_expr_vec > 4.0] = 4.0
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_", marker_gene, "_expression.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(marker_gene_expr_vec, umap_loadings, marker_gene)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

marker_gene = "CD14"
marker_gene_index <- which(gene_names==marker_gene)
marker_gene_expr_vec <- expr[, marker_gene_index]
marker_gene_expr_vec[marker_gene_expr_vec > 4.0] = 4.0
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_", marker_gene, "_expression.pdf")
umap_scatter_cd14 <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(marker_gene_expr_vec, umap_loadings, marker_gene)
ggsave(umap_scatter_cd14, file=output_file, width=7.2, height=6.0, units="in")


marker_gene = "RGS1"
marker_gene_index <- which(gene_names==marker_gene)
marker_gene_expr_vec <- expr[, marker_gene_index]
marker_gene_expr_vec[marker_gene_expr_vec > 4.0] = 4.0
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_", marker_gene, "_expression.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(marker_gene_expr_vec, umap_loadings, marker_gene)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

marker_gene = "NKG7"
marker_gene_index <- which(gene_names==marker_gene)
marker_gene_expr_vec <- expr[, marker_gene_index]
marker_gene_expr_vec[marker_gene_expr_vec > 4.0] = 4.0
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_", marker_gene, "_expression.pdf")
umap_scatter_nkg7 <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(marker_gene_expr_vec, umap_loadings, marker_gene)
ggsave(umap_scatter_nkg7, file=output_file, width=7.2, height=6.0, units="in")

marker_gene = "PRF1"
marker_gene_index <- which(gene_names==marker_gene)
marker_gene_expr_vec <- expr[, marker_gene_index]
marker_gene_expr_vec[marker_gene_expr_vec > 4.0] = 4.0
output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_", marker_gene, "_expression.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_variable(marker_gene_expr_vec, umap_loadings, marker_gene)
ggsave(umap_scatter, file=output_file, width=7.2, height=6.0, units="in")

output_file <- paste0(visualization_dir, "umap_loading_scatter_colored_by_multiple_marker_genes_expression.pdf")
merged <- plot_grid(umap_scatter_cd14, umap_scatter_nkg7, umap_scatter_cd8b, umap_scatter_bank1, ncol=2, labels = c('A', 'B', 'C', 'D'))
ggsave(merged, file=output_file, width=7.2, height=8.0, units="in")


