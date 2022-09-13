args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(cowplot)
library(umap)
library(ggplot2)
library(RColorBrewer)
options(bitmapType = 'cairo', device = 'pdf')


make_loading_scatter_plot_for_fixed_dimensions <- function(loadings, tissues, factor_1, factor_2, tissue_colors) {
	unique_tissues = unique(tissues)
	colors <- c()
	for (tissue_iter in 1:length(unique_tissues)) {
		tiss <- unique_tissues[tissue_iter]
		hex = tissue_colors$tissue_color_hex[tissue_colors$tissue_id == tiss]
		colors <- c(colors, paste0("#",hex))

	}


	df <- data.frame(loading_1=loadings[,factor_1], loading_2=loadings[,factor_2], tissue=factor(tissues))
	plotter <- ggplot(df) + 
	           geom_point( aes(x=loading_1, y=loading_2, color=tissue),size=.01) +
	           gtex_v8_figure_theme() + 
	           scale_color_manual(values=colors) +
	           labs(x=paste0("Loading ", factor_1), y = paste0("Loading ", factor_2), color="") +
	           theme(legend.position="none")
	return(plotter)
}

make_loading_scatter_plot_for_fixed_dimensions_legend <- function(loadings, tissues, factor_1, factor_2, tissue_colors) {
	unique_tissues = unique(tissues)
	colors <- c()
	for (tissue_iter in 1:length(unique_tissues)) {
		tiss <- unique_tissues[tissue_iter]
		hex = tissue_colors$tissue_color_hex[tissue_colors$tissue_id == tiss]
		colors <- c(colors, paste0("#",hex))

	}

	df <- data.frame(loading_1=loadings[,factor_1], loading_2=loadings[,factor_2], tissue=factor(tissues))
	plotter <- ggplot(df) + 
	           geom_point( aes(x=loading_1, y=loading_2, color=tissue),size=.01) +
	           gtex_v8_figure_theme() + 
	           scale_color_manual(values=colors) +
	           labs(x=paste0("Loading ", factor_1), y = paste0("Loading ", factor_2), color="") +
	           guides(colour=guide_legend(nrow=2,byrow=TRUE, override.aes = list(size=2))) +
	           theme(legend.position="bottom")

	return(get_legend(plotter))
}

make_loading_scatter_plot <- function(tissues, sample_covariate_file, tissue_colors, loading_file) {
	#tissues <- read.table(tissue_file, header=FALSE)
	loadings <- read.table(loading_file, header=FALSE)
	covariates <- read.table(sample_covariate_file, header=TRUE, sep="\t")
	race_factor = factor(as.character(as.numeric(covariates$race==3)))
	num_samples <- dim(loadings)[1]
	num_factors <- dim(loadings)[2]
	if (num_factors == 2) {
	unique_tissues = unique(tissues)
	colors <- c()
	for (tissue_iter in 1:length(unique_tissues)) {
		tiss <- unique_tissues[tissue_iter]
		hex = tissue_colors$tissue_color_hex[tissue_colors$tissue_id == tiss]
		colors <- c(colors, paste0("#",hex))

	}
	df <- data.frame(loading_1=loadings$V1, loading_2=loadings$V2, tissue=factor(tissues))
	plotter <- ggplot(df) + 
	           geom_point( aes(x=loading_1, y=loading_2, color=tissue),size=.005) +
	           gtex_v8_figure_theme() + 
	           scale_color_manual(values=colors) + 
	           labs(x="Loading 1", y = "Loading 2", color="")
	} else if (num_factors == 3) {
		plot_arr <- c()
		list_counter <- 1
		for (factor_1 in 1:num_factors) {
			for (factor_2 in 1:num_factors) {
				temp_plot <- make_loading_scatter_plot_for_fixed_dimensions(loadings, tissues, factor_1, factor_2, tissue_colors)
				if (factor_1 == factor_2) {
					plot_arr[[list_counter]] = NULL
				} else {
					plot_arr[[list_counter]] <- temp_plot
				}
				
				list_counter <- list_counter + 1

			}
		}

		legend <- make_loading_scatter_plot_for_fixed_dimensions_legend(loadings, tissues, factor_1, factor_2, tissue_colors)
		
		plotter <- plot_grid(plot_grid(plot_arr[[1]], plot_arr[[2]], plot_arr[[3]], plot_arr[[4]], plot_arr[[5]], plot_arr[[6]], plot_arr[[7]], plot_arr[[8]], ncol=num_factors), legend, ncol=1,rel_heights=c(1,.15))
	} else if (num_factors == 4) {
		plot_arr <- c()
		list_counter <- 1
		for (factor_1 in 1:num_factors) {
			for (factor_2 in 1:num_factors) {
				temp_plot <- make_loading_scatter_plot_for_fixed_dimensions(loadings, tissues, factor_1, factor_2, tissue_colors)
				if (factor_1 == factor_2) {
					plot_arr[[list_counter]] = NULL
				} else {
					plot_arr[[list_counter]] <- temp_plot
				}
				
				list_counter <- list_counter + 1

			}
		}

		legend <- make_loading_scatter_plot_for_fixed_dimensions_legend(loadings, tissues, factor_1, factor_2, tissue_colors)
		
		plotter <- plot_grid(plot_grid(plot_arr[[2]], plot_arr[[3]], plot_arr[[4]], plot_arr[[6]], plot_arr[[7]], plot_arr[[8]], NULL, plot_arr[[11]], plot_arr[[12]], ncol=num_factors-1), legend, ncol=1,rel_heights=c(1,.15))

	} else if (num_factors == 5) {
		plot_arr <- c()
		list_counter <- 1
		for (factor_1 in 1:num_factors) {
			for (factor_2 in 1:num_factors) {
				temp_plot <- make_loading_scatter_plot_for_fixed_dimensions(loadings, tissues, factor_1, factor_2, tissue_colors)
				if (factor_1 == factor_2) {
					plot_arr[[list_counter]] = NULL
				} else {
					plot_arr[[list_counter]] <- temp_plot
				}
				
				list_counter <- list_counter + 1

			}
		}

		legend <- make_loading_scatter_plot_for_fixed_dimensions_legend(loadings, tissues, factor_1, factor_2, tissue_colors)
		
		plotter <- plot_grid(plot_grid(plot_arr[[1]], plot_arr[[2]], plot_arr[[3]], plot_arr[[4]], plot_arr[[5]], plot_arr[[6]], plot_arr[[7]], plot_arr[[8]], plot_arr[[9]], plot_arr[[10]], plot_arr[[11]], plot_arr[[12]], plot_arr[[13]], plot_arr[[14]], plot_arr[[15]], plot_arr[[16]],plot_arr[[17]], plot_arr[[18]], plot_arr[[19]], plot_arr[[20]], plot_arr[[21]], plot_arr[[22]], plot_arr[[23]], plot_arr[[24]], ncol=num_factors), legend, ncol=1,rel_heights=c(1,.15))

	}

	return(plotter)
}

make_umap_loading_scatter_plot <- function(tissues, tissue_colors, sample_covariate_file, loading_file, umap_loadings) {
	covariates <- read.table(sample_covariate_file, header=TRUE, sep="\t")
	loadings <- read.table(loading_file, header=FALSE)

	unique_tissues = unique(tissues)
	colors <- c()
	for (tissue_iter in 1:length(unique_tissues)) {
		tiss <- unique_tissues[tissue_iter]
		hex = tissue_colors$tissue_color_hex[tissue_colors$tissue_id == tiss]
		colors <- c(colors, paste0("#",hex))

	}
	# 0.034993038124563364, -0.02809246824633377, 0.06527282210541974
	# umap_loadings = umap(loadings)$layout

	df <- data.frame(loading_1=umap_loadings[,1], loading_2=umap_loadings[,2], tissue=factor(tissues), race=factor(as.character(as.numeric(covariates$race==3))))
	plotter <- ggplot(df) + 
	           geom_point(aes(x=loading_1, y=loading_2, color=tissue, shape=race), size=1) +
	           scale_color_manual(values=colors) + 
	           scale_shape_manual(values = c(15,3)) +
	           gtex_v8_figure_theme() + 
	           guides(colour = guide_legend(override.aes = list(size=2))) +
	           labs(x="UMAP 1", y = "UMAP 2", color="", shape="Race") + 
	           guides(colour=guide_legend(nrow=5,byrow=TRUE, override.aes = list(size=2))) +
	           theme(legend.position="bottom") + 
	           theme(legend.text = element_text(size=8), legend.title = element_text(size=8))
	return(plotter)
}


make_umap_loading_scatter_plot_colored_by_real_valued_vector <- function(umap_loadings, real_valued_covariate, covariate_name) {
	df <- data.frame(loading_1=umap_loadings[,1], loading_2=umap_loadings[,2], covariate=real_valued_covariate)
	plotter <- ggplot(df) + 
	           geom_point(aes(x=loading_1, y=loading_2, color=covariate), size=1) +
	           gtex_v8_figure_theme() + 
	           labs(x="UMAP 1", y = "UMAP 2", color=covariate_name) + 
	           theme(legend.position="bottom")
	return(plotter)
}

make_loading_boxplot_plot_by_race <- function(sample_covariate_file, loadings) {
	#tissues <- read.table(tissue_file, header=FALSE)
	#loadings <- read.table(loading_file, header=FALSE)
	#df <- data.frame(loading_1=loadings$V1, loading_2=loadings$V2, tissue=factor(tissues$V1))
	covariates <- read.table(sample_covariate_file, header=TRUE, sep="\t")

	loading_vec <- c()
	race_vec <- c()
	factor_vec <- c()

	num_factors <- dim(loadings)[2]

	for (factor_number in 1:num_factors) {
		loading_vec <- c(loading_vec, loadings[,factor_number])
		#race_vec <- c(race_vec, as.character(as.numeric(covariates$race==3)))
		race_vec <- c(race_vec, covariates$race)
		factor_vec <- c(factor_vec, rep(as.character(factor_number), length(covariates$race)))
	}


	df <- data.frame(loading=loading_vec, race=factor(race_vec), latent_factor=factor(factor_vec, levels=as.character(1:num_factors)))

	df = df[(df$race==3) | (df$race==2),]

	ancestry <- c()
	num_samples = length(df$race)
	for (sample_num in 1:num_samples) {
		if (df$race[sample_num] == 3) {
			ancestry <- c(ancestry, "European Ancestry")
		} else {
			ancestry <- c(ancestry, "African Ancestry")
		}
	}
	df$ancestry = factor(ancestry)


	boxplot <- ggplot(df, aes(x=latent_factor, y=loading, fill=ancestry)) + geom_boxplot(outlier.size = .1) +
				gtex_v8_figure_theme() + 
	        	labs(x="eQTL Factorization latent context", y = "Measurement loading", fill="") +
	        	theme(legend.position="bottom")

	return(boxplot)
}

make_loading_boxplot_plot_by_cm_proportion <- function(sample_covariate_file, loading_file) {
	#tissues <- read.table(tissue_file, header=FALSE)
	loadings <- read.table(loading_file, header=FALSE)
	#df <- data.frame(loading_1=loadings$V1, loading_2=loadings$V2, tissue=factor(tissues$V1))
	covariates <- read.table(sample_covariate_file, header=TRUE, sep="\t")

	loading_vec <- c()
	race_vec <- c()
	factor_vec <- c()

	num_factors <- dim(loadings)[2]

	for (factor_number in 1:num_factors) {
		loading_vec <- c(loading_vec, loadings[,factor_number])
		race_vec <- c(race_vec, as.character(as.numeric(covariates$cm_cell_type_composition < .5)))
		factor_vec <- c(factor_vec, rep(as.character(factor_number), length(covariates$cm_cell_type_composition)))
	}


	df <- data.frame(loading=loading_vec, race=factor(race_vec), latent_factor=factor(factor_vec, levels=as.character(1:num_factors)))

	boxplot <- ggplot(df, aes(x=latent_factor, y=loading, fill=race)) + geom_boxplot(outlier.size = .1) +
				gtex_v8_figure_theme() + 
	        	labs(x="Latent factor", y = "Sample loading", fill="Known Cell type composition") +
	        	theme(legend.position="bottom")

	return(boxplot)
}

make_loading_boxplot_plot_by_sex <- function(sample_covariate_file, loadings) {
	#tissues <- read.table(tissue_file, header=FALSE)
	#loadings <- read.table(loading_file, header=FALSE)
	#df <- data.frame(loading_1=loadings$V1, loading_2=loadings$V2, tissue=factor(tissues$V1))
	covariates <- read.table(sample_covariate_file, header=TRUE, sep="\t")

	loading_vec <- c()
	sex_vec <- c()
	factor_vec <- c()

	num_factors <- dim(loadings)[2]

	for (factor_number in 1:num_factors) {
		loading_vec <- c(loading_vec, loadings[,factor_number])
		sex_vec <- c(sex_vec, as.character(as.numeric(covariates$sex==2)))
		factor_vec <- c(factor_vec, rep(as.character(factor_number), length(covariates$sex)))
	}


	df <- data.frame(loading=loading_vec, sex=factor(sex_vec), latent_factor=factor(factor_vec, levels=as.character(1:num_factors)))

	boxplot <- ggplot(df, aes(x=latent_factor, y=loading, fill=sex)) + geom_boxplot(outlier.size = .1) +
				gtex_v8_figure_theme() + 
	        	labs(x="Latent factor", y = "Sample loading", fill="Known sex") +
	        	theme(legend.position="bottom")

	return(boxplot)
}

make_loading_boxplot_plot_by_cohort <- function(sample_covariate_file, loadings) {
	#tissues <- read.table(tissue_file, header=FALSE)
	#loadings <- read.table(loading_file, header=FALSE)
	#df <- data.frame(loading_1=loadings$V1, loading_2=loadings$V2, tissue=factor(tissues$V1))
	covariates <- read.table(sample_covariate_file, header=TRUE, sep="\t")

	loading_vec <- c()
	cohort_vec <- c()
	factor_vec <- c()

	num_factors <- dim(loadings)[2]

	for (factor_number in 1:num_factors) {
		loading_vec <- c(loading_vec, loadings[,factor_number])
		cohort_vec <- c(cohort_vec, as.character(covariates$cohort=="Postmortem"))
		factor_vec <- c(factor_vec, rep(as.character(factor_number), length(covariates$cohort)))
	}


	df <- data.frame(loading=loading_vec, cohort=factor(cohort_vec), latent_factor=factor(factor_vec, levels=as.character(1:num_factors)))

	boxplot <- ggplot(df, aes(x=latent_factor, y=loading, fill=cohort)) + geom_boxplot(outlier.size = .1) +
				gtex_v8_figure_theme() + 
	        	labs(x="Latent factor", y = "Sample loading", fill="Known cohort") +
	        	theme(legend.position="bottom")

	return(boxplot)
}

make_loading_boxplot_plot_by_age <- function(sample_covariate_file, loading_file) {
	#tissues <- read.table(tissue_file, header=FALSE)
	loadings <- read.table(loading_file, header=FALSE)
	#df <- data.frame(loading_1=loadings$V1, loading_2=loadings$V2, tissue=factor(tissues$V1))
	covariates <- read.table(sample_covariate_file, header=TRUE, sep="\t")

	loading_vec <- c()
	race_vec <- c()
	factor_vec <- c()

	num_factors <- dim(loadings)[2]

	for (factor_number in 1:num_factors) {
		loading_vec <- c(loading_vec, loadings[,factor_number])
		age_vec <- c(race_vec, covariates$age)
		print(factor_number)
		print(cor(loadings[,factor_number],age_vec))
		factor_vec <- c(factor_vec, rep(as.character(factor_number), length(covariates$age)))
	}

	if (FALSE) {
	df <- data.frame(loading=loading_vec, sex=factor(sex_vec), latent_factor=factor(factor_vec, levels=as.character(1:num_factors)))

	boxplot <- ggplot(df, aes(x=latent_factor, y=loading, fill=sex)) + geom_boxplot(outlier.size = .1) +
				gtex_v8_figure_theme() + 
	        	labs(x="Latent factor", y = "Sample loading", fill="Known sex") +
	        	theme(legend.position="bottom")

	return(boxplot)
	}
}

make_loading_boxplot_plot_by_tissue <- function(tissues,tissue_colors, loadings) {
	#tissues <- read.table(tissue_file, header=FALSE)
	#loadings <- read.table(loading_file, header=FALSE)
	#df <- data.frame(loading_1=loadings$V1, loading_2=loadings$V2, tissue=factor(tissues$V1))

	loading_vec <- c()
	tissue_vec <- c()
	factor_vec <- c()

	num_factors <- dim(loadings)[2]

	for (factor_number in 1:num_factors) {
		loading_vec <- c(loading_vec, loadings[,factor_number])
		tissue_vec <- c(tissue_vec, as.character(tissues))
		factor_vec <- c(factor_vec, rep(as.character(factor_number), length(tissues)))
	}


	df <- data.frame(loading=loading_vec, tissue=factor(tissue_vec), latent_factor=factor(factor_vec, levels=as.character(1:num_factors)))
	unique_tissues = unique(df$tissue)
	colors <- c()
	for (tissue_iter in 1:length(unique_tissues)) {
		tiss <- unique_tissues[tissue_iter]
		hex = tissue_colors$tissue_color_hex[tissue_colors$tissue_id == tiss]
		colors <- c(colors, paste0("#",hex))

	}

	boxplot <- ggplot(df, aes(x=latent_factor, y=loading, fill=tissue)) + geom_boxplot(outlier.size = .001) +
				gtex_v8_figure_theme() + 
				scale_fill_manual(values=colors) + 
	        	labs(x="eQTL Factorization latent context", y = "Measurement loading", fill="") +
	        	theme(legend.position="bottom") +
	           	guides(fill=guide_legend(nrow=4,byrow=TRUE, override.aes = list(size=.1))) + 
	           	theme(legend.text=element_text(size=9))

	return(boxplot)
}


make_loading_boxplot_plot_by_tissue_non_color <- function(tissues, loading_file) {
	#tissues <- read.table(tissue_file, header=FALSE)
	loadings <- read.table(loading_file, header=FALSE)
	#df <- data.frame(loading_1=loadings$V1, loading_2=loadings$V2, tissue=factor(tissues$V1))

	loading_vec <- c()
	tissue_vec <- c()
	factor_vec <- c()

	num_factors <- dim(loadings)[2]

	for (factor_number in 1:num_factors) {
		loading_vec <- c(loading_vec, loadings[,factor_number])
		tissue_vec <- c(tissue_vec, as.character(tissues))
		factor_vec <- c(factor_vec, rep(as.character(factor_number), length(tissues)))
	}


	df <- data.frame(loading=loading_vec, tissue=factor(tissue_vec), latent_factor=factor(factor_vec, levels=as.character(1:num_factors)))


	boxplot <- ggplot(df, aes(x=latent_factor, y=loading, fill=tissue)) + geom_boxplot(outlier.size = .001) +
				gtex_v8_figure_theme() + 
	        	labs(x="eQTL Factorization latent context", y = "Measurement loading", fill="") +
	        	theme(legend.position="bottom") +
	           	guides(fill=guide_legend(nrow=4,byrow=TRUE, override.aes = list(size=.1))) + 
	           	theme(legend.text=element_text(size=9))

	return(boxplot)
}


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}

get_tissue_names <- function(sample_file_name) {
	aa <- read.table(sample_file_name)
	vecy = as.character(aa$V1)
	tissues <- c()
	for (iter in 1:length(vecy)) {
		tissue <- strsplit(vecy[iter],":")[[1]][2]
		tissues <- c(tissues, tissue)
	}
	return(tissues)
}


get_indi_names <- function(sample_file_name) {
	aa <- read.table(sample_file_name)
	vecy = as.character(aa$V1)
	tissues <- c()
	for (iter in 1:length(vecy)) {
		tissue <- strsplit(vecy[iter],":")[[1]][1]
		tissues <- c(tissues, tissue)
	}
	return(tissues)
}

 make_residual_clustering_heatmap <- function(lm_residual_file, tissue_names) {
 	# Get matrix of dimension num_genesXnum_indi
 	residual_mat <- abs(t(read.table(lm_residual_file)))
 	#residual_mat <- residual_mat[1:100, 1:200]
 	residual_mat[residual_mat > 3] = 3.0
 	num_genes <- dim(residual_mat)[1]
 	num_indi <- dim(residual_mat)[2]
	colnames(residual_mat) = paste0("indi_", 1:num_indi)
	rownames(residual_mat) = paste0("gene_", 1:num_genes)

	ord <- hclust( dist(scale(residual_mat), method = "euclidean"), method = "ward.D" )$order

	melted_mat <- melt(residual_mat)
    colnames(melted_mat) <- c("gene", "sample", "error")

    #melted_mat$gene <- factor(paste0("gene_", 1:num_genes), levels = paste0("gene_", 1:num_genes)[ord])
    melted_mat$gene <- factor(melted_mat$gene, levels = paste0("gene_", 1:num_genes)[ord])

    melted_mat$sample <- factor(melted_mat$sample, levels = paste0("indi_", 1:num_indi))


    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=gene, y=sample)) + geom_tile(aes(fill=error)) + 
		gtex_v8_figure_theme() +
   		labs(y="RNA-seq Sample", x="Variant-gene", fill="Absolute lm residual") +
   		theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
   		theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
   		scale_fill_gradient(low="pink",high="blue") +
   		scale_x_discrete(expand=c(0,0)) +
 		scale_y_discrete(expand=c(0,0))
	return(heatmap)
 }

  make_abs_expr_clustering_heatmap <- function(expr_file, tissue_names) {
  	print(expr_file)
 	# Get matrix of dimension num_genesXnum_indi
 	residual_mat <- abs(t(t(read.table(expr_file))))
 	residual_mat[residual_mat > 3] = 3.0
 	#residual_mat[residual_mat < -3] = -3.0
 	num_genes <- dim(residual_mat)[1]
 	num_indi <- dim(residual_mat)[2]


	colnames(residual_mat) = paste0("indi_", 1:num_indi)
	rownames(residual_mat) = paste0("gene_", 1:num_genes)


	ord <- hclust( dist(scale(residual_mat), method = "euclidean"), method = "ward.D" )$order
	melted_mat <- melt(residual_mat)
    colnames(melted_mat) <- c("gene", "sample", "error")

    #melted_mat$gene <- factor(paste0("gene_", 1:num_genes), levels = paste0("gene_", 1:num_genes)[ord])
    melted_mat$gene <- factor(melted_mat$gene, levels = paste0("gene_", 1:num_genes)[ord])

    melted_mat$sample <- factor(melted_mat$sample, levels = paste0("indi_", 1:num_indi))


    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=gene, y=sample)) + geom_tile(aes(fill=error)) + 
		gtex_v8_figure_theme() +
   		labs(y="RNA-seq Sample", x="Variant-gene", fill="Abs Expr") +
   		theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
   		theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
   		scale_fill_gradient(low="pink",high="blue") +
   		scale_x_discrete(expand=c(0,0)) +
 		scale_y_discrete(expand=c(0,0))
	return(heatmap)
 }

 make_expr_clustering_heatmap <- function(expr_file, tissue_names) {
  	print(expr_file)
 	# Get matrix of dimension num_genesXnum_indi
 	residual_mat <- t(t(read.table(expr_file)))
 	residual_mat[residual_mat > 3] = 3.0
 	residual_mat[residual_mat < -3] = -3.0
 	num_genes <- dim(residual_mat)[1]
 	num_indi <- dim(residual_mat)[2]


	colnames(residual_mat) = paste0("indi_", 1:num_indi)
	rownames(residual_mat) = paste0("gene_", 1:num_genes)


	ord <- hclust( dist(scale(residual_mat), method = "euclidean"), method = "ward.D" )$order
	melted_mat <- melt(residual_mat)
    colnames(melted_mat) <- c("gene", "sample", "error")

    #melted_mat$gene <- factor(paste0("gene_", 1:num_genes), levels = paste0("gene_", 1:num_genes)[ord])
    melted_mat$gene <- factor(melted_mat$gene, levels = paste0("gene_", 1:num_genes)[ord])

    melted_mat$sample <- factor(melted_mat$sample, levels = paste0("indi_", 1:num_indi))


    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=gene, y=sample)) + geom_tile(aes(fill=error)) + 
		gtex_v8_figure_theme() +
   		labs(y="RNA-seq Sample", x="Variant-gene", fill="Expr") +
   		theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
   		theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
   		scale_fill_gradient(low="pink",high="blue") +
   		scale_x_discrete(expand=c(0,0)) +
 		scale_y_discrete(expand=c(0,0))
	return(heatmap)
 }

elbo_line_plot_for_various_random_initializations <- function(eqtl_results_dir, model_stem, num_seeds, hot_start) {
	elbo_arr <- c()
	seed_id_arr <- c()
	iter_arr <- c()
	for (seed_number in 0:(num_seeds-1)) {
		seed_elbo_file <- paste0(eqtl_results_dir, model_stem , seed_number, "_seed_elbo.txt")
		temp_data <- read.table(seed_elbo_file, header=FALSE)
		elbo = temp_data[,1]
		if (hot_start == "True") {
			elbo = elbo[5:length(elbo)]
		}
		elbo_arr <- c(elbo_arr, elbo)
		seed_id_arr <- c(seed_id_arr, rep(paste0("seed ", seed_number), length(elbo)))
		iter_arr <- c(iter_arr, 1:length(elbo))
	}
	df <- data.frame(elbo=elbo_arr, iterations=iter_arr, seed=factor(seed_id_arr))


	p<-ggplot(df, aes(x=iterations, y=elbo, group=seed)) +
  		geom_line(aes(color=seed))+
		gtex_v8_figure_theme() +
		labs(x="Variational inference iteration", y="ELBO", color="")
	return(p)
}

make_absolute_effect_size_boxplot <- function(effect_size_file, tissue_colors, factor_file, num_tests) {
	effect_sizes <- read.table(effect_size_file, header=TRUE)
	all_factors <- read.table(factor_file, header=FALSE)
	num_tissues = dim(effect_sizes)[2]
	num_factors <- dim(all_factors)[1]

	effect_size_arr <- c()
	factor_num_arr <- c()
	tissue_name_arr <- c()

	unique_tissues = unique(colnames(effect_sizes))
	colors <- c()
	for (tissue_iter in 1:length(unique_tissues)) {
		tiss <- unique_tissues[tissue_iter]
		hex = tissue_colors$tissue_color_hex[tissue_colors$tissue_id == tiss]
		colors <- c(colors, paste0("#",hex))

	}

	for (factor_num in 1:num_factors) {
		# Find indices of top n factors
		sorted_factor <- sort(abs(all_factors[factor_num,]), decreasing=TRUE)
		min_val = as.numeric(sorted_factor[num_tests])
		indices = abs(all_factors[factor_num,]) >= min_val
	
		# Put info into data frame
		for (tissue_num in 1:num_tissues) {
			tissue_name <- colnames(effect_sizes)[tissue_num]
			tissue_effect_sizes <- effect_sizes[indices, tissue_num]
			effect_size_arr <- c(effect_size_arr, abs(tissue_effect_sizes))
			tissue_name_arr <- c(tissue_name_arr, rep(tissue_name, length(tissue_effect_sizes)))
			factor_num_arr <- c(factor_num_arr, rep(factor_num, length(tissue_effect_sizes)))
		}

	}
	df <- data.frame(effect_size=effect_size_arr, tissue=factor(tissue_name_arr), latent_factor=factor(factor_num_arr))
	boxplot <- ggplot(df, aes(x=latent_factor, y=effect_size, fill=tissue)) + geom_violin() +
				gtex_v8_figure_theme() + 
				scale_fill_manual(values=colors) + 
	        	labs(x="Latent factor", y = "Absolute effect size", fill="Known tissue") +
	        	guides(colour=guide_legend(nrow=4,byrow=TRUE, override.aes = list(size=2))) +
	        	theme(legend.position="bottom")

	return(boxplot)

}


#Make heatmap showing correlation of loadings between models
make_loading_correlation_heatmap <- function(model1_loading_file, model2_loading_file, x_axis_label, y_axis_label) {

	# Load in loading matrices
	model1_loadings <- read.table(model1_loading_file, header=FALSE)
	model2_loadings <- read.table(model2_loading_file, header=FALSE)

	num_dim_x = dim(model1_loadings)[2]
	num_dim_y = dim(model2_loadings)[2]
	
	corr_matrix = matrix(0, num_dim_x, num_dim_y)

	for (x_index in 1:num_dim_x) {
		for (y_index in 1:num_dim_y) {
			corr_matrix[x_index, y_index] <- abs(cor(model1_loadings[,x_index], model2_loadings[, y_index]))
		}
	}
	
	melted_mat <- melt(corr_matrix)
    colnames(melted_mat) <- c("model1", "model2", "correlation")


    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=model1, y=model2)) + geom_tile(aes(fill=correlation)) + 
		gtex_v8_figure_theme() +
   		labs(y=y_axis_label, x=x_axis_label, fill="Absolute\nPearson correlation") +
   		#theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
   		#theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
   		scale_fill_gradient(low="white",high="blue") +
   		scale_x_continuous(breaks=1:num_dim_x, labels=1:num_dim_x) +
   		scale_y_continuous(breaks=1:num_dim_y, labels=1:num_dim_y)
	return(heatmap)

}


make_ancestry_scatter_plot_in_one_tissue <- function(tissue, ancestry_specific_factor, ancestry_specific_eqtl_effect_size_file) {
	ancestry_specific_eqtl_effect_sizes <- read.table(ancestry_specific_eqtl_effect_size_file, header=TRUE)

	ancestry_specific_eqtl_effect_sizes$aa_effect_size[ancestry_specific_eqtl_effect_sizes$aa_effect_size == "NaN"] = NA
	ancestry_specific_eqtl_effect_sizes$ea_effect_size[ancestry_specific_eqtl_effect_sizes$ea_effect_size == "NaN"] = NA

	# Put into cmpact data frame
	df <- data.frame(factor=abs(ancestry_specific_factor), delta_effect_size=abs(ancestry_specific_eqtl_effect_sizes$ea_effect_size-ancestry_specific_eqtl_effect_sizes$aa_effect_size))

	test_scores <- cor.test(df$factor, df$delta_effect_size, use="complete.obs")
	plotter <- ggplot(df, aes(x=delta_effect_size, y=factor)) + 
	           geom_point(size=.01, colour="slateblue3") +
	           gtex_v8_figure_theme() + 
	           labs(x="Abs(European Ancestry - African Ancestry cis-eQTL effect size)", y ="Abs(Factor weight)", title=paste0(tissue, " / spearman rho: ", signif(test_scores$estimate, 4), " / p-value: ", signif(test_scores$p.value, 4))) 
	return(plotter)

}
make_ancestry_scatter_plot <- function(factor_file, tissue_names, ancestry_specific_eqtl_file_root) {
	# Extract ancestry specific factor
	factors <- read.table(factor_file, header=FALSE)
	ancestry_specific_factor <- as.numeric(factors[4,])
	# make scatter plot in 1 tissue
	tissue_1_ancestry_scatter <- make_ancestry_scatter_plot_in_one_tissue(tissue_names[1], ancestry_specific_factor, paste0(ancestry_specific_eqtl_file_root, tissue_names[1], "_effect_sizes.txt"))

	tissue_2_ancestry_scatter <- make_ancestry_scatter_plot_in_one_tissue(tissue_names[2], ancestry_specific_factor, paste0(ancestry_specific_eqtl_file_root, tissue_names[2], "_effect_sizes.txt"))

	tissue_3_ancestry_scatter <- make_ancestry_scatter_plot_in_one_tissue(tissue_names[3], ancestry_specific_factor, paste0(ancestry_specific_eqtl_file_root, tissue_names[3], "_effect_sizes.txt"))

	tissue_4_ancestry_scatter <- make_ancestry_scatter_plot_in_one_tissue(tissue_names[4], ancestry_specific_factor, paste0(ancestry_specific_eqtl_file_root, tissue_names[4], "_effect_sizes.txt"))

	combined_ancestry_scatter <- plot_grid(tissue_1_ancestry_scatter, tissue_2_ancestry_scatter, tissue_3_ancestry_scatter, tissue_4_ancestry_scatter, ncol=1)

	return(combined_ancestry_scatter)
}

make_ancestry_missingness_boxplot <- function(factor_file, tissue_names, ancestry_specific_eqtl_file_root) {
	# Extract ancestry specific factor
	factors <- read.table(factor_file, header=FALSE)
	ancestry_specific_factor <- as.numeric(factors[4,])
	# Extract ancestry specific effect sizes
	ancestry_specific_eqtl_effect_sizes_1 <- read.table(paste0(ancestry_specific_eqtl_file_root, tissue_names[1], "_effect_sizes.txt"), header=TRUE)
	ancestry_specific_eqtl_effect_sizes_2 <- read.table(paste0(ancestry_specific_eqtl_file_root, tissue_names[2], "_effect_sizes.txt"), header=TRUE)
	ancestry_specific_eqtl_effect_sizes_3 <- read.table(paste0(ancestry_specific_eqtl_file_root, tissue_names[3], "_effect_sizes.txt"), header=TRUE)
	ancestry_specific_eqtl_effect_sizes_4 <- read.table(paste0(ancestry_specific_eqtl_file_root, tissue_names[4], "_effect_sizes.txt"), header=TRUE)

	delta_effect_size_1 <- ancestry_specific_eqtl_effect_sizes_1$ea_effect_size - ancestry_specific_eqtl_effect_sizes_1$aa_effect_size
	delta_effect_size_2 <- ancestry_specific_eqtl_effect_sizes_2$ea_effect_size - ancestry_specific_eqtl_effect_sizes_2$aa_effect_size
	delta_effect_size_3 <- ancestry_specific_eqtl_effect_sizes_3$ea_effect_size - ancestry_specific_eqtl_effect_sizes_3$aa_effect_size
	delta_effect_size_4 <- ancestry_specific_eqtl_effect_sizes_4$ea_effect_size - ancestry_specific_eqtl_effect_sizes_4$aa_effect_size

	num_tests <- length(ancestry_specific_factor)
	median_delta_effect_sizes <- c()

	for (test_index in 1:num_tests) {
		med_delta_effect_size <- median(c(delta_effect_size_1[test_index], delta_effect_size_2[test_index], delta_effect_size_3[test_index], delta_effect_size_4[test_index]))
		median_delta_effect_sizes <- c(median_delta_effect_sizes, med_delta_effect_size)
	}
	
	df <- data.frame(factor=abs(ancestry_specific_factor), missingness=is.na(median_delta_effect_sizes))

	#test_scores <- cor.test(df$factor, df$delta_effect_size, use="complete.obs")
	plotter <- ggplot(df, aes(x=factor(missingness), y=factor)) + 
				geom_boxplot(outlier.size = .1) +
	           gtex_v8_figure_theme() + 
	           labs(x="Missingness of eQTLs in an ancestry", y ="Abs(Factor weight)") 
	return(plotter)
}

make_median_ancestry_scatter_plot <- function(factor_file, tissue_names, ancestry_specific_eqtl_file_root) {
	# Extract ancestry specific factor
	factors <- read.table(factor_file, header=FALSE)
	ancestry_specific_factor <- as.numeric(factors[4,])
	# Extract ancestry specific effect sizes
	ancestry_specific_eqtl_effect_sizes_1 <- read.table(paste0(ancestry_specific_eqtl_file_root, tissue_names[1], "_effect_sizes.txt"), header=TRUE)
	ancestry_specific_eqtl_effect_sizes_2 <- read.table(paste0(ancestry_specific_eqtl_file_root, tissue_names[2], "_effect_sizes.txt"), header=TRUE)
	ancestry_specific_eqtl_effect_sizes_3 <- read.table(paste0(ancestry_specific_eqtl_file_root, tissue_names[3], "_effect_sizes.txt"), header=TRUE)
	ancestry_specific_eqtl_effect_sizes_4 <- read.table(paste0(ancestry_specific_eqtl_file_root, tissue_names[4], "_effect_sizes.txt"), header=TRUE)

	delta_effect_size_1 <- ancestry_specific_eqtl_effect_sizes_1$ea_effect_size - ancestry_specific_eqtl_effect_sizes_1$aa_effect_size
	delta_effect_size_2 <- ancestry_specific_eqtl_effect_sizes_2$ea_effect_size - ancestry_specific_eqtl_effect_sizes_2$aa_effect_size
	delta_effect_size_3 <- ancestry_specific_eqtl_effect_sizes_3$ea_effect_size - ancestry_specific_eqtl_effect_sizes_3$aa_effect_size
	delta_effect_size_4 <- ancestry_specific_eqtl_effect_sizes_4$ea_effect_size - ancestry_specific_eqtl_effect_sizes_4$aa_effect_size

	num_tests <- length(ancestry_specific_factor)
	median_delta_effect_sizes <- c()

	for (test_index in 1:num_tests) {
		med_delta_effect_size <- median(c(delta_effect_size_1[test_index], delta_effect_size_2[test_index], delta_effect_size_3[test_index], delta_effect_size_4[test_index]))
		median_delta_effect_sizes <- c(median_delta_effect_sizes, med_delta_effect_size)
	}
	
	df <- data.frame(factor=abs(ancestry_specific_factor), delta_effect_size=abs(median_delta_effect_sizes))

	test_scores <- cor.test(df$factor, df$delta_effect_size, use="complete.obs")
	plotter <- ggplot(df, aes(x=delta_effect_size, y=factor)) + 
	           geom_point(size=.01, colour="slateblue3") +
	           gtex_v8_figure_theme() + 
	           labs(x="Abs(Median(European Ancestry - African Ancestry cis-eQTL effect size))", y ="Abs(Factor weight)", title=paste0("spearman rho: ", signif(test_scores$estimate, 4), " / p-value: ", signif(test_scores$p.value, 4))) 
	return(plotter)
}


make_cell_type_loading_scatter <- function(loading, cell_type_enrichment, loading_name, cell_type_name, tissues, tissue_colors) {
	df <- data.frame(loading=loading, cell_type_enrichment=cell_type_enrichment, tissue=factor(tissues))

	unique_tissues = unique(df$tissue)
	colors <- c()
	for (tissue_iter in 1:length(unique_tissues)) {
		tiss <- unique_tissues[tissue_iter]
		hex = tissue_colors$tissue_color_hex[tissue_colors$tissue_id == tiss]
		colors <- c(colors, paste0("#",hex))

	}
	plotter <- ggplot(df, aes(x=loading, y=cell_type_enrichment, color=tissue)) + 
	           geom_point() +
	           gtex_v8_figure_theme() + 
	           scale_color_manual(values=colors) +
	           theme(legend.position="none") +	
	           labs(x=loading_name, y=cell_type_name)

	return(plotter)

}

make_cell_type_loadings_scatters <- function(loadings, cell_type_enrichment, cell_type_name, tissues, tissue_colors) {
	num_factors <- dim(loadings)[2]

	plot_arr <- list()
	for (loading_num in 1:num_factors) {
		p <- make_cell_type_loading_scatter(loadings[,loading_num], cell_type_enrichment, paste0("Loading ", loading_num), cell_type_name, tissues, tissue_colors) 
		plot_arr[[loading_num]] <- p
	}
	merged = plot_grid(plotlist=plot_arr, ncol=2)

	return(merged)

}

make_cell_type_loading_scatter_no_gtex_color <- function(loading, cell_type_enrichment, loading_name, cell_type_name, cov, cov_name) {
	df <- data.frame(loading=loading, cell_type_enrichment=cell_type_enrichment, cov=factor(cov))


	plotter <- ggplot(df, aes(x=loading, y=cell_type_enrichment, color=cov)) + 
	           geom_point() +
	           gtex_v8_figure_theme() + 
	           labs(x=loading_name, y=cell_type_name, color=cov_name)

	return(plotter)

}


make_cell_type_loadings_scatters_no_gtex_color <- function(loadings, cell_type_enrichment, cell_type_name, cov, cov_name) {
	num_factors <- dim(loadings)[2]

	plot_arr <- list()
	for (loading_num in 1:num_factors) {
		p <- make_cell_type_loading_scatter_no_gtex_color(loadings[,loading_num], cell_type_enrichment, paste0("Loading ", loading_num), cell_type_name, cov, cov_name) 
		plot_arr[[loading_num]] <- p
	}
	merged = plot_grid(plotlist=plot_arr, ncol=2)

	return(merged)

}

make_factor_distribution_histograms <- function(tissue_10_factor_file) {
	factors <- as.matrix(read.table(tissue_10_factor_file, header=FALSE))

	factor_min <- floor(min(factors)) 
	factor_max <- ceiling(max(factors))

	num_factors <- dim(factors)[1]

	plot_arr <- list()
	factor_vals_vec <- c()
	factor_names_vec <- c()
	factor_num <- 1
	num_tests <- dim(factors)[2]
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

make_cell_type_loadings_scatter_for_samples_from_specified_tissue <- function(loading, cell_type_enrichment, cell_type_name, tissue_name, loading_name, tissues, valid_tissues, tissue_colors) {
	df <- data.frame(loading=loading, cell_type_enrichment=cell_type_enrichment, tissue=factor(tissues))
	df2 = df[as.character(df$tissue) %in% as.character(valid_tissues),]

	unique_tissues = unique(df2$tissue)
	colors <- c()
	for (tissue_iter in 1:length(unique_tissues)) {
		tiss <- unique_tissues[tissue_iter]
		hex = tissue_colors$tissue_color_hex[tissue_colors$tissue_id == tiss]
		colors <- c(colors, paste0("#",hex))

	}

	corry <- cor.test(df2$loading, df2$cell_type_enrichment)
	plotter <- ggplot(df2, aes(x=loading, y=cell_type_enrichment, color=tissue)) + 
			   scale_color_manual(values=colors) +
	           geom_point() +
	           gtex_v8_figure_theme() + 
	           labs(x=loading_name, y=cell_type_name, color="", title=paste0("Pearson correlation pvalue: ", signif(corry$p.value, digits=4))) + 
	           guides(color=guide_legend(nrow=4,byrow=TRUE))

	return(plotter)
}

make_cell_type_loadings_scatter_for_samples_from_specified_tissue_colored_by_categorical_variable <- function(loading, cell_type_enrichment, cell_type_name, tissue_name, loading_name, tissues, valid_tissues, color_variable, color_variable_name) {
	df <- data.frame(loading=loading, cell_type_enrichment=cell_type_enrichment, tissue=factor(tissues), cat_variable=factor(color_variable))
	df2 = df[as.character(df$tissue) %in% as.character(valid_tissues),]


	corry <- cor.test(df2$loading, df2$cell_type_enrichment)
	plotter <- ggplot(df2, aes(x=loading, y=cell_type_enrichment, color=cat_variable)) + 
	           geom_point() +
	           gtex_v8_figure_theme() + 
	           labs(x=loading_name, y=cell_type_name, title=paste0(tissue_name, " / correlation pvalue: ", corry$p.value), color=color_variable_name)

	return(plotter)

}


make_cell_type_loadings_scatter_colored_by_specified_tissues_and_categorical_variables <- function(loading, cell_type_enrichment, cell_type_name, tissue_name, loading_name, tissues, valid_tissues, color_variable, color_variable_name) {
	df <- data.frame(loading=loading, cell_type_enrichment=cell_type_enrichment, tissue=factor(tissues), cat_variable=factor(color_variable))
	df$tissue_subset = factor(as.character(df$tissue) %in% as.character(valid_tissues), levels=c(TRUE, FALSE))



	#corry <- cor.test(df2$loading, df2$cell_type_enrichment)
	plotter <- ggplot(df, aes(x=loading, y=cell_type_enrichment, color=tissue_subset, shape=cat_variable)) + 
	           geom_point(size=1.2, alpha=.5) +
	           gtex_v8_figure_theme() + 
	           labs(x=loading_name, y=cell_type_name, color=paste0(tissue_name), shape=color_variable_name)

	return(plotter)

}

make_cell_type_loadings_scatter_colored_by_tissues_and_categorical_variables <- function(loading, cell_type_enrichment, cell_type_name, loading_name, tissues, tissue_colors, valid_tissues, categorical_variable, categorical_variable_name) {
	df <- data.frame(loading=loading, cell_type_enrichment=cell_type_enrichment, tissue=factor(tissues), cat_variable=factor(categorical_variable))

	unique_tissues = unique(df$tissue)
	colors <- c()
	for (tissue_iter in 1:length(unique_tissues)) {
		tiss <- unique_tissues[tissue_iter]
		hex = tissue_colors$tissue_color_hex[tissue_colors$tissue_id == tiss]
		colors <- c(colors, paste0("#",hex))

	}

	plotter <- ggplot(df, aes(x=loading, y=cell_type_enrichment, color=tissue, shape=cat_variable)) + 
			   geom_point() + 
	           gtex_v8_figure_theme() + 
	           scale_color_manual(values=colors) +
	           labs(x=loading_name, y=cell_type_name, color="", shape=categorical_variable_name)

	return(plotter)

}



make_loadings_loadings_scatter_colored_by_cell_type_for_samples_from_specified_tissue <- function(loadings1, loadings2, cell_type_enrichment, cell_type_name, tissue_name, loading_1_name, loading_2_name, tissues, valid_tissues) {
	df <- data.frame(loading1=loadings1, loading2=loadings2, cell_type_enrichment=cell_type_enrichment, tissue=factor(tissues))
	df2 = df[as.character(df$tissue) %in% as.character(valid_tissues),]
	plotter <- ggplot(df2, aes(x=loading1, y=loading2, color=cell_type_enrichment)) + 
	           geom_point() +
	           gtex_v8_figure_theme() + 
	           labs(x=loading_1_name, y=loading_2_name, color=cell_type_name, title=paste0(tissue_name))

	return(plotter)
}

loading_by_genotype_pc1_colored_by_race <- function(loading,race, genotype_pc) {
	df <- data.frame(loading=loading, race=race, genotype_pc=genotype_pc)
	plotter <- ggplot(df, aes(x=loading, y=genotype_pc, color=race)) + 
	           geom_point() +
	           gtex_v8_figure_theme() + 
	           labs(x="loading 1", y="genotype PC1", color="Known Ancestry")
}




make_loading_boxplot_plot_by_batch_for_single_factor <- function(technical_covariate_file, loading_file, loading_num) {
	loadings <- read.table(loading_file, header=FALSE)
	covs <- read.table(technical_covariate_file, header=TRUE, sep="\t")

	df <- data.frame(loading=loadings[, loading_num], batch=covs$SMGEBTCH)

	plotter <- ggplot(df, aes(x=batch, y=loading)) + 
		       geom_boxplot(outlier.size = .1) +
	           gtex_v8_figure_theme() + 
	           labs(x="Batch", y =paste0("Loading ", loading_num)) +
	           theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
	return(plotter) 


}


make_loading_boxplot_plot_by_individual_for_single_factor <- function(indi_names, loading_file, loading_num) {
	loadings <- read.table(loading_file, header=FALSE)
	#covs <- read.table(technical_covariate_file, header=TRUE, sep="\t")

	df <- data.frame(loading=loadings[, loading_num], batch=indi_names)

	plotter <- ggplot(df, aes(x=batch, y=loading)) + 
		       geom_boxplot(outlier.size = .1) +
	           gtex_v8_figure_theme() + 
	           labs(x="Donor ID", y =paste0("Loading ", loading_num)) +
	           theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
	return(plotter) 


}

explore_relationship_between_surveyed_covariates_and_eqtl_factors <- function(surveyed_covariate_file, loadings) {
	#loadings <- read.table(loading_file, header=FALSE)
	covs <- read.table(surveyed_covariate_file, header=TRUE)

	print(summary(covs))
	covs <- covs[,2:(dim(covs)[2])]
	factor_num <- 2
	num_covs = dim(covs)[2]
	vec <- c()
	for (cov_num in 1:num_covs) {
		cov_vec <- covs[, cov_num]
		loading_vec <- loadings[,factor_num]

		valid_samples <- cov_vec == 0.0 | cov_vec == 1
		lin_model <- lm(loading_vec[valid_samples] ~ cov_vec[valid_samples])
		print(paste0(cov_num, "    ", summary(lin_model)$adj.r.squared))
		vec <- c(vec, summary(lin_model)$adj.r.squared)
	}
	print(max(vec))
	#print(summary(covs))
}

explore_relationship_between_sample_covariates_and_eqtl_factors <- function(sample_covariate_file, loadings) {
	#loadings <- read.table(loading_file, header=FALSE)
	covs <- read.table(sample_covariate_file, header=TRUE, sep="\t")
	covs <- covs[,2:(dim(covs)[2])]
	factor_num <- 2
	num_covs = dim(covs)[2]
	loading_vec <- loadings[, factor_num]
	print(colnames(covs))

	for (cov_num in 2:num_covs) {
		cov_vec <- covs[, cov_num]
		loading_vec <- loadings[,factor_num]

		lin_model <- lm(loading_vec ~ cov_vec)
		print(paste0(cov_num, "    ", summary(lin_model)$adj.r.squared, "    ", summary(lin_model)$r.squared))
	}
}

explore_relationship_between_muscle_age_and_eqtl_factors <- function(sample_covariate_file, loadings) {
	#loadings <- read.table(loading_file, header=FALSE)
	covs <- read.table(sample_covariate_file, header=TRUE, sep="\t")
	covs <- covs[,2:(dim(covs)[2])]

	factor_num <- 1
	loading_vec <- loadings[, factor_num]

	muscle_tissue_binary = 1.0*(as.character(covs$tissue_type) == "Muscle_Skeletal")


	lin_model <- lm(loading_vec ~ muscle_tissue_binary + covs$age + covs$age:muscle_tissue_binary)

	print(summary(lin_model))
}


explore_relationship_between_technical_covariates_and_eqtl_factors <- function(surveyed_covariate_file, loadings, indi_names) {
	#loadings <- read.table(loading_file, header=FALSE)
	covs <- read.table(surveyed_covariate_file, header=TRUE, sep="\t")
	covs <- covs[,2:(dim(covs)[2])]
	print(summary(covs))
	#factor_num <- 8
	factor_num <- 2
	num_covs = dim(covs)[2]
	loading_vec <- loadings[, factor_num]
	indi_names <- factor(indi_names)

	lin_model <- lm(loading_vec ~ indi_names)
	print(summary(lin_model)$adj.r.squared)
	print(summary(lin_model)$r.squared)

	for (cov_num in 1:num_covs) {
		cov_vec <- covs[, cov_num]
		loading_vec <- loadings[,factor_num]

		lin_model <- lm(loading_vec ~ cov_vec)
		print(paste0(cov_num, "    ", summary(lin_model)$adj.r.squared, "    ", summary(lin_model)$r.squared))
	}
	#print(summary(covs))
}

#Make heatmap showing correlation of loadings between models
make_eqtl_factor_expression_pc_correlation_heatmap <- function(loadings, expression_pc_file, x_axis_label, y_axis_label) {

	# Load in loading matrices
	pcs <- read.table(expression_pc_file, header=TRUE)
	# Remove sample name column from pcs
	pcs <- pcs[,2:(dim(pcs)[2])]

	num_dim_x = dim(loadings)[2]
	num_dim_y = dim(pcs)[2]
	
	corr_matrix = matrix(0, num_dim_x, num_dim_y)

	for (x_index in 1:num_dim_x) {
		for (y_index in 1:num_dim_y) {
			if (var(loadings[,x_index]) == 0.0) {
				corr_matrix[x_index, y_index] <- NA
			} else {
				corr_matrix[x_index, y_index] <- abs(cor(loadings[,x_index], pcs[, y_index]))
			}
		}
	}
	

	colnames(corr_matrix) = colnames(pcs)
	rownames(corr_matrix) = paste0(1:num_dim_x)


	melted_mat <- melt(corr_matrix)
    colnames(melted_mat) <- c("model1", "model2", "correlation")


    melted_mat$model2 <- factor(melted_mat$model2, levels=colnames(pcs))
    melted_mat$model1 <- factor(melted_mat$model1, levels=paste0(1:num_dim_x))


    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=model1, y=model2)) + geom_tile(aes(fill=correlation)) + 
		gtex_v8_figure_theme() +
   		labs(y=y_axis_label, x=x_axis_label, fill="Absolute\nPearson correlation") +
   		scale_fill_gradient(low="white",high="blue") 
	return(heatmap)

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
                labs(x = "Factor number", y = "Fraction of QTL variance explained") + 
                gtex_v8_figure_theme() 

    return(line_plot)
}
make_loading_boxplot_for_one_factor_by_race_and_tissue_type <- function(race,tissue_names, loadings, factor_number) {
  df <- data.frame(loading=loadings, race=factor(race), tissue=tissue_names)

  df = df[(df$race==3) | (df$race==2),]

  	ancestry <- c()
	num_samples = length(df$race)
	for (sample_num in 1:num_samples) {
		if (df$race[sample_num] == 3) {
			ancestry <- c(ancestry, "European Ancestry")
		} else {
			ancestry <- c(ancestry, "African Ancestry")
		}
	}
	df$ancestry = factor(ancestry)

  boxplot <- ggplot(df, aes(x=tissue, y=loading, fill=ancestry)) + geom_boxplot(outlier.size = .00001) +
        gtex_v8_figure_theme() + 
            labs(x="", y = paste0("Sample loading (", factor_number,")"), fill="") +
            theme(legend.position="top") +
            guides(colour = guide_legend(override.aes = list(size=2))) 
     
   return(boxplot)

}


make_loading_boxplot_plot_with_row_for_every_factor_by_race_and_tissue_type<- function(race, tissue_names, loadings) {
  loading_vec <- c()
  covariate_vec <- c()
  num_factors <- dim(loadings)[2]

  plot_arr <- list()

  legend <- get_legend(make_loading_boxplot_for_one_factor_by_race_and_tissue_type(race, tissue_names, loadings[,num_factors], num_factors))

  plot_arr[[1]] = legend

  for (factor_num in 1:(num_factors-1)) {
    factor_boxplot <- make_loading_boxplot_for_one_factor_by_race_and_tissue_type(race, tissue_names, loadings[,factor_num], factor_num) + theme(axis.text.x=element_blank()) + theme(legend.position="none")
    plot_arr[[(1+factor_num)]] <- factor_boxplot
  }

  plot_arr[[(1+num_factors)]] <- make_loading_boxplot_for_one_factor_by_race_and_tissue_type(race, tissue_names, loadings[,num_factors], num_factors) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position="none")



  merged = plot_grid(plotlist=plot_arr, ncol=1, rel_heights=c(.6,1,1,1,1,1,2.9))


  return(merged)
}

make_loading_by_cell_type_scatter <- function(loading_vec, cell_type_vec, loading_num, cell_type_name) {
  df <- data.frame(loading=loading_vec, cell_type_prop=cell_type_vec)
  lm1 = lm(loading~cell_type_prop, data=df)
  r_squared = summary(lm1)$adj.r.squared

plotter <- ggplot(df, aes(x=loading, y=cell_type_vec)) + 
           geom_point() +
           gtex_v8_figure_theme() + 
           labs(x=paste0("SURGE Loading ", loading_num), y=cell_type_name, title=paste0("Adjusted r squard: ", r_squared))
    return(plotter)

}
make_loading_by_lm_cell_type_scatter <- function(loading_vec, loading_number, covariates) {
	lm1 = lm(loading_vec~Adipocytes+Epithelial_cells+Hepatocytes+Keratinocytes+Myocytes+Neurons+Neutrophils, data=covariates) #Create the linear regression
	p <- predict(lm1, covariates)
  df <- data.frame(loading=loading_vec, cell_type_prop=p)

  r_squared = summary(lm1)$adj.r.squared

plotter <- ggplot(df, aes(x=loading, y=cell_type_prop)) + 
           geom_point() +
           gtex_v8_figure_theme() + 
           labs(x=paste0("Loading ", loading_number), y="cell-type proportion lm", title=paste0("Adjusted r squard: ", r_squared))
    return(plotter)

}

make_loading_cell_type_pvalue_plot <- function(loading_vec, loading_number, covariates) {
	lm1 = lm(loading_vec~Adipocytes+Epithelial_cells+Hepatocytes+Keratinocytes+Myocytes+Neurons+Neutrophils, data=covariates) #Create the linear regression
 }

make_cell_type_enrichment_pvalue_plot <- function(loading_vec, covariates, loading_number) {
	#print(scale(covariates))
	new_df <- data.frame(Adipocytes=scale(covariates$Adipocytes), Epithelial_cells=scale(covariates$Epithelial_cells), Hepatocytes=scale(covariates$Hepatocytes), Keratinocytes=scale(covariates$Keratinocytes), Myocytes=scale(covariates$Myocytes), Neurons=scale(covariates$Neurons), Neutrophils=scale(covariates$Neutrophils))

	lm1 = lm(loading_vec~Adipocytes+Epithelial_cells+Hepatocytes+Keratinocytes+Myocytes+Neurons+Neutrophils, data=new_df) #Create the linear regression
	
	coefs = coef(summary(lm1))[2:(dim(coef(summary(lm1)))[1]), ]
	ub = coefs[,1] + 2.96*coefs[,2]
	lb = coefs[,1] - 2.96*coefs[,2]
	df <- data.frame(term=row.names(coefs), estimate=coefs[,1], conf.low=lb, conf.high=ub)
	p <- ggplot(df, aes(term, estimate,color=term))+
  		geom_point()+
  		geom_pointrange(aes(ymin = conf.low, ymax = conf.high))+
  		gtex_v8_figure_theme() +
  		geom_hline(yintercept=0) +
  		labs(x="", y="Effect size") +
  		theme(legend.position="none") +
  		theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  	return(p)
}
observed_cell_type_proportion_stacked_bar_plot <- function(loading_vec, loading_number, covariates, num_bins) {
	minny <- min(loading_vec)
	maxy <- max(loading_vec)
	bin_vec = rep("NULL", length(loading_vec))
	cutoff_vec <- seq(minny, maxy, length.out=(num_bins+1))
	ordered_bins <- c()
	for (bin_num in 1:num_bins) {
		lower = cutoff_vec[bin_num]
		upper = cutoff_vec[(bin_num + 1)]
		if (bin_num < num_bins) {
			bin_indices = loading_vec >= lower & loading_vec < upper 
			bin_vec[bin_indices] = rep(paste0("bin", bin_num), sum(bin_indices))
		} else {
			bin_indices = loading_vec >= lower & loading_vec <= upper
			bin_vec[bin_indices] = rep(paste0("bin", bin_num), sum(bin_indices))
		}
		ordered_bins <- c(ordered_bins, paste0(bin_num))


		# debug
		#temp_indices = bin_vec == paste0("bin", bin_num)
	}


	cell_type_vec <- c()
	final_bin_vec <- c()
	count_vec <- c()
	ordered_cell_types <- c("Adipocytes", "Epithelial_cells", "Hepatocytes", "Keratinocytes", "Myocytes", "Neurons", "Neutrophils")

	for (cell_type_iter in 1:length(ordered_cell_types)) {
		for (bin_num in 1:num_bins) {
			cell_type <- ordered_cell_types[cell_type_iter]

			bin_indices <- bin_vec == paste0("bin", bin_num)
			temp_cell_type_specific_counts <- covariates[[cell_type]]

			counts = sum(temp_cell_type_specific_counts[bin_indices])

			cell_type_vec <- c(cell_type_vec, cell_type)
			final_bin_vec <- c(final_bin_vec, paste0(bin_num))
			count_vec <- c(count_vec, counts)

		}
	}

	df <- data.frame(cell_type=factor(cell_type_vec, ordered_cell_types), bin=factor(final_bin_vec, levels=ordered_bins), counts=count_vec)
	# Stacked + percent
	p <- ggplot(df, aes(fill=cell_type, y=counts, x=bin)) + 
    	geom_bar(position="fill", stat="identity") +
    	gtex_v8_figure_theme() +
    	labs(fill="", y="Cell type proportions", x=paste0("SURGE Loading ", loading_number, " bin"))
    return(p)
}


loading_by_various_cell_types_scatters <- function(loading_vec, loading_number, covariates) {
	adipo_scatter <- make_loading_by_cell_type_scatter(loading_vec, covariates$Adipocytes, loading_number, "Adipocytes")
	epi_scatter <- make_loading_by_cell_type_scatter(loading_vec, covariates$Epithelial_cells, loading_number, "Epithelial")
	hep_scatter <- make_loading_by_cell_type_scatter(loading_vec, covariates$Hepatocytes, loading_number, "Hepatocytes")
	ker_scatter <- make_loading_by_cell_type_scatter(loading_vec, covariates$Keratinocytes, loading_number, "Keratinocytes")
	myo_scatter <- make_loading_by_cell_type_scatter(loading_vec, covariates$Myocytes, loading_number, "Myocytes")
	nuer_scatter <- make_loading_by_cell_type_scatter(loading_vec, covariates$Neurons, loading_number, "Neurons")
	neut_scatter <- make_loading_by_cell_type_scatter(loading_vec, covariates$Neutrophils, loading_number, "Neutrophils")
	pred_scatter <- make_loading_by_lm_cell_type_scatter(loading_vec, loading_number, covariates)

	#pred_scatter <- make_loading_cell_type_pvalue_plot(loading_vec, loading_number, covariates)

	merged = plot_grid(adipo_scatter, epi_scatter, hep_scatter, ker_scatter, myo_scatter, nuer_scatter, neut_scatter, pred_scatter, ncol=2)

	return(merged)
}

sample_non_significant_hits <- function(pvalues, fraction_kept=.01,fraction_sampled=.001) {
    index <- floor(length(pvalues)*fraction_kept)
    to_keep <- pvalues[1:index]
    to_filter <- pvalues[(index+1):length(pvalues)]
    filtered <- sort(sample(to_filter,floor(length(to_filter)*fraction_sampled)))
    return(c(to_keep,filtered))
}

make_interaction_eqtl_qq_plot <- function(interaction_eqtl_dir, tissue_stem) {
	surge_qtl_file <- paste0(interaction_eqtl_dir, tissue_stem, "_surge_results_k_init_10_seed_1_warmup_5_ratio_variance_std_True_perm_False_interaction_eqtl_results_latent_factor_1_merged.txt")
	surge_qtl_pvalues <- read.table(surge_qtl_file, header=FALSE)$V5
	sorted_surge_qtl_pvalues <- sample_non_significant_hits(sort(surge_qtl_pvalues))

	sorted_null_qtl_pvalues <- sample_non_significant_hits(sort(runif(length(surge_qtl_pvalues))))
	


	df <- data.frame(surge_pvalue=-log10(sorted_surge_qtl_pvalues+1e-40), cell_type_pvalue=-log10(sorted_null_qtl_pvalues+1e-40))

    # PLOT!
    max_val <-max(max(-log10(sorted_surge_qtl_pvalues + 1e-40)), max(-log10(sorted_null_qtl_pvalues + 1e-40)))
    #PLOT!
    scatter <- ggplot(df, aes(x = cell_type_pvalue, y = surge_pvalue)) + geom_point()
    scatter <- scatter + theme(text = element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
    scatter <- scatter + labs(colour="",x = expression(log[10]("Expected eQTL p-value")), y = expression(log[10]("SURGE interaction eQTL p-value")))
    scatter <- scatter + geom_abline()
    scatter <- scatter + theme(legend.position="bottom")
    scatter <- scatter + scale_x_continuous(limits = c(-.1, max_val + .1), breaks = round(seq(0, max_val, by = 5),1)) + scale_y_continuous(limits = c(-.1,max_val+.1), breaks = round(seq(0, max_val, by = 5),1))
    scatter <- scatter + theme(plot.title=element_text(size=8, face="plain"), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8))
    return(scatter)

}


make_interaction_eqtl_cell_type_comparison_qq_plot <- function(interaction_eqtl_dir, tissue_stem) {
	surge_qtl_file <- paste0(interaction_eqtl_dir, tissue_stem, "_surge_results_k_init_10_seed_1_warmup_5_ratio_variance_std_True_perm_False_interaction_eqtl_results_latent_factor_1_merged.txt")
	surge_qtl_pvalues <- read.table(surge_qtl_file, header=FALSE)$V5
	sorted_surge_qtl_pvalues <- sample_non_significant_hits(sort(surge_qtl_pvalues))

	cell_types <- c("Adipocytes", "Epithelial_cells", "Hepatocytes", "Keratinocytes", "Myocytes", "Neurons", "Neutrophils")
	
	surge_pvalue_vec <- c()
	cell_type_pvalue_vec <- c()
	cell_type_vec <- c()
	for (cell_type_iter in 1:length(cell_types)) {
		cell_type <- cell_types[cell_type_iter]
		cell_type_qtl_file <- paste0(interaction_eqtl_dir, "xcell_interaction_tissues_subset_colon_transverse_", cell_type, "_interaction_eqtl_results_merged_pvalue_only.txt")
		cell_type_qtl_pvalues <- read.table(cell_type_qtl_file, header=FALSE)$V1
		sorted_cell_type_qtl_pvalues <- sort(cell_type_qtl_pvalues)



		surge_pvalue_vec <- c(surge_pvalue_vec, sorted_surge_qtl_pvalues)
		temper = sample_non_significant_hits(sorted_cell_type_qtl_pvalues)
		cell_type_pvalue_vec <- c(cell_type_pvalue_vec, temper)
		cell_type_vec <- c(cell_type_vec, rep(cell_type, length(sorted_surge_qtl_pvalues)))
	}
	df <- data.frame(surge_pvalue=-log10(surge_pvalue_vec+1e-100), cell_type_pvalue=-log10(cell_type_pvalue_vec+1e-100), cell_type=factor(cell_type_vec, levels=cell_types))

    # PLOT!
    max_val <-max(max(-log10(surge_pvalue_vec + 1e-100)), max(-log10(cell_type_pvalue_vec + 1e-100)))
    #PLOT!
    scatter <- ggplot(df, aes(x = cell_type_pvalue, y = surge_pvalue, colour = cell_type)) + geom_point()
    scatter <- scatter + theme(text = element_text(size=14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
    scatter <- scatter + labs(colour="",x = expression(log[10]("Cell type interaction eQTL p-value")), y = expression(log[10]("SURGE interaction eQTL p-value")))
    scatter <- scatter + geom_abline()
    scatter <- scatter + theme(legend.position="bottom")
    scatter <- scatter + scale_x_continuous(limits = c(-.1, max_val + .1), breaks = round(seq(0, max_val, by = 5),1)) + scale_y_continuous(limits = c(-.1,max_val+.1), breaks = round(seq(0, max_val, by = 5),1))
    scatter <- scatter + theme(plot.title=element_text(size=8, face="plain"), text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8))
    return(scatter)

}




processed_data_dir <- args[1]
eqtl_results_dir <- args[2]
visualization_dir <- args[3]
tissue_colors_file <- args[4]
interaction_eqtl_dir <- args[5]


# Read in tissue colors and names
tissue_colors = read.table(tissue_colors_file, header = T, stringsAsFactors = F, sep = "\t")

tissue_colors$tissue_id = tissue_colors$tissue_site_detail_id


############################
# Load in files
############################
options(warn=1)
stem <- "tissues_subset_colon_transverse"
tissue_10_file <- paste0(processed_data_dir, stem, "_outliers_removed_sample_names.txt")
tissue_10_sample_covariate_file <- paste0(processed_data_dir, stem, "_outliers_removed_sample_covariates.txt")
tissue_10_surveyed_covariate_file <- paste0(processed_data_dir, stem, "_outliers_removed_sample_surveyed_covariates.txt")
tissue_10_technical_covariate_file <- paste0(processed_data_dir, stem, "_outliers_removed_sample_technical_covariates.txt")
tissue_10_expression_pcs_file <- paste0(processed_data_dir, stem, "_covariates.txt")
#tissue_10_residual_expression_pcs_file <- paste0(processed_data_dir, "tissues_subset_10_v2_residual_expression_covariates_120_pc.txt")

tissue_10_names <- get_tissue_names(tissue_10_file)
tissue_10_indi_names <- get_indi_names(tissue_10_file)


############################
# Model Specification
############################
tissue_10_model_stem <- paste0(stem, "_surge_results_k_init_10_seed_1_warmup_5_ratio_variance_std_True_permute_False_2000_tests_var_param_1e-3_")
tissue_10_loading_file <- paste0(eqtl_results_dir, tissue_10_model_stem, "U_S.txt")
tissue_10_factor_file <- paste0(eqtl_results_dir, tissue_10_model_stem, "V.txt")

pve_file <- paste0(eqtl_results_dir, tissue_10_model_stem, "factor_pve.txt")
pve <- as.numeric(read.table(pve_file, header=FALSE, sep="\t")$V1)
ordering <- order(pve, decreasing=TRUE)
#ordering <- ordering[1:8]



loadings <- read.table(tissue_10_loading_file, header=FALSE)
loadings <- loadings[, ordering]
ordered_pve <- pve[ordering]

############################
# Start making plots!!
############################

# Load in loading matrices
pcs <- read.table(tissue_10_expression_pcs_file, header=TRUE)
# Remove sample name column from pcs
pcs <- pcs[,2:(dim(pcs)[2])]


covariates <- read.table(tissue_10_sample_covariate_file, header=TRUE, sep="\t")

#######################################
# Make interaction eqtl QQ plot
#######################################
output_file <- paste0(visualization_dir, tissue_10_model_stem, "interaction_eqtl_qq_plot.pdf")
qq_plot <- make_interaction_eqtl_qq_plot(interaction_eqtl_dir, stem)
ggsave(qq_plot, file=output_file, width=7.2, height=5.5, units="in")



output_file <- paste0(visualization_dir, tissue_10_model_stem, "interaction_eqtl_cell_type_comparison_qq_plot.pdf")
qq_plot <- make_interaction_eqtl_cell_type_comparison_qq_plot(interaction_eqtl_dir, stem)
ggsave(qq_plot, file=output_file, width=7.2, height=5.5, units="in")


#######################################
# PVE plot showing fraction of eqtl variance explained through each factor
#######################################
output_file <- paste0(visualization_dir, tissue_10_model_stem, "fraction_of_eqtl_variance_explained_lineplot.pdf")
pve_plot <- make_pc_variance_explained_line_plot(ordered_pve)
ggsave(pve_plot, file=output_file, width=7.2, height=5.5, units="in")

##################
# Stacked barplot of precicted cell type proportions along latent factor
loading_number <- 1
num_bins <- 10
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_", loading_number, "_observed_", num_bins, "_bins_cell_type_proportions_stacked_barplot.pdf")
boxplots <- observed_cell_type_proportion_stacked_bar_plot(loadings[,loading_number], loading_number, covariates, num_bins)
ggsave(boxplots, file=output_file, width=7.2, height=4.0, units="in")


##################
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_1_by_cell_type_enrichment_pvalues.pdf")
cell_type_pvalues_plot <- make_cell_type_enrichment_pvalue_plot(loadings[,1], covariates, 1)
ggsave(cell_type_pvalues_plot, file=output_file, width=7.2, height=5.5, units="in")



joint_plot <- plot_grid(boxplots, cell_type_pvalues_plot, labels=c("A","B"), ncol=1)
output_file <- paste0(visualization_dir, tissue_10_model_stem, "supplentary_cell_type_composition.pdf")
ggsave(joint_plot, file=output_file, width=7.2, height=6.0, units="in")


































##################
# OLD
##################





if (FALSE) {

#####################
# Make heatmap correlation eqtl factorization loadings with expression pcs
output_file <- paste0(visualization_dir, tissue_10_model_stem, "expression_pc_eqtl_factorization_loading_correlation_heatmap.pdf")
#heatmap <- make_eqtl_factor_expression_pc_correlation_heatmap(loadings, tissue_10_expression_pcs_file, "eQTL Factorization Loadings", "Expression PC Loadings")
#ggsave(heatmap, file=output_file, width=7.2, height=8.5, units="in")


######################
# Make box plot for each Race, showing loading distributions
output_file <- paste0(visualization_dir, tissue_10_model_stem, "race_colored_loading_boxplot.pdf")
boxplot_race <- make_loading_boxplot_plot_by_race(tissue_10_sample_covariate_file, loadings)
ggsave(boxplot_race, file=output_file, width=7.2, height=5.5, units="in")


######################
# Make box plot for each Cohort, showing loading distributions
output_file <- paste0(visualization_dir, tissue_10_model_stem, "cohort_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_cohort(tissue_10_sample_covariate_file, loadings)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")

######################
# Make box plot for each Sex, showing loading distributions
output_file <- paste0(visualization_dir, tissue_10_model_stem, "sex_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_sex(tissue_10_sample_covariate_file, loadings)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")

######################
# Make box plot for each tissue, showing loading distributions
output_file <- paste0(visualization_dir, tissue_10_model_stem, "tissue_colored_loading_boxplot.pdf")
boxplot_tissue <- make_loading_boxplot_plot_by_tissue(tissue_10_names, tissue_colors, loadings)
ggsave(boxplot_tissue, file=output_file, width=7.2, height=5.5, units="in")

######################
# Make cowplot merged boxplot
output_file <- paste0(visualization_dir, tissue_10_model_stem, "tissue_and_ancestry_colored_loading_boxplot.pdf")
tissue_legend <- get_legend(boxplot_tissue)
race_legend <- get_legend(boxplot_race)
merged <- plot_grid(boxplot_tissue + theme(legend.position="none"), boxplot_race + theme(legend.position="none"), tissue_legend, race_legend, ncol=2, rel_heights=c(1,.5), rel_widths=c(1, .8), labels = c('a', 'b', '', ''))
ggsave(merged, file=output_file, width=13.2, height=5.5, units="in")

##################
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_by_genotype_pc1_colored_by_race.pdf")
boxplots <- loading_by_genotype_pc1_colored_by_race(loadings[,1],factor(covariates$race), pcs$genotype_PC0)
ggsave(boxplots, file=output_file, width=7.2, height=5.5, units="in")


##################
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_1_by_epithelial_cell_enrichment.pdf")
epi_scatter <- make_loading_by_cell_type_scatter(loadings[,1], covariates$Epithelial_cells, 1, "Epithelial cell enrichment")
ggsave(epi_scatter, file=output_file, width=7.2, height=5.5, units="in")




}

##################
# Stacked barplot of precicted cell type proportions along latent factor
loading_number <- 1
num_bins <- 10
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_", loading_number, "_observed_", num_bins, "_bins_cell_type_proportions_stacked_barplot.pdf")
boxplots <- observed_cell_type_proportion_stacked_bar_plot(loadings[,loading_number], loading_number, covariates, num_bins)
ggsave(boxplots, file=output_file, width=7.2, height=4.0, units="in")

if (FALSE) {

loading_number <- 2
num_bins <- 10
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_", loading_number, "_observed_", num_bins, "_bins_cell_type_proportions_stacked_barplot.pdf")
boxplots <- observed_cell_type_proportion_stacked_bar_plot(loadings[,loading_number], loading_number, covariates, num_bins)
ggsave(boxplots, file=output_file, width=7.2, height=4.0, units="in")



loading_number <- 2
num_bins <- 10
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_", loading_number, "_observed_", num_bins, "_bins_cell_type_proportions_stacked_barplot.pdf")
if (var(loadings[,loading_number]) != 0.0) {
	boxplots <- observed_cell_type_proportion_stacked_bar_plot(loadings[,loading_number], loading_number, covariates, num_bins)
	ggsave(boxplots, file=output_file, width=7.2, height=4.0, units="in")
}

loading_number <- 3
num_bins <- 10
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_", loading_number, "_observed_", num_bins, "_bins_cell_type_proportions_stacked_barplot.pdf")
if (var(loadings[,loading_number]) != 0.0) {
	boxplots <- observed_cell_type_proportion_stacked_bar_plot(loadings[,loading_number], loading_number, covariates, num_bins)
	ggsave(boxplots, file=output_file, width=7.2, height=4.0, units="in")
}

##################
# Make scatter plot showing correlation between loadings and cell type enrichments
loading_number <- 1
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_", loading_number, "_by_various_cell_types_scatter_plots.pdf")
boxplots <- loading_by_various_cell_types_scatters(loadings[,loading_number], loading_number, covariates)
ggsave(boxplots, file=output_file, width=7.2, height=8.5, units="in")

##################
# Make scatter plot showing correlation between loadings and cell type enrichments
loading_number <- 2
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_", loading_number, "_by_various_cell_types_scatter_plots.pdf")
boxplots <- loading_by_various_cell_types_scatters(loadings[,loading_number], loading_number, covariates)
ggsave(boxplots, file=output_file, width=7.2, height=8.5, units="in")

##################
# Make scatter plot showing correlation between loadings and cell type enrichments
loading_number <- 3
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_", loading_number, "_by_various_cell_types_scatter_plots.pdf")
boxplots <- loading_by_various_cell_types_scatters(loadings[,loading_number], loading_number, covariates)
ggsave(boxplots, file=output_file, width=7.2, height=8.5, units="in")


##################
# Make scatter plot showing correlation between loadings and cell type enrichments
output_file <- paste0(visualization_dir, tissue_10_model_stem, "line_per_factor_loadings_boxplot_seperated_by_ancestry_and_tissue.pdf")
boxplots <- make_loading_boxplot_plot_with_row_for_every_factor_by_race_and_tissue_type(covariates$race, tissue_10_names, loadings[,1:6])
ggsave(boxplots, file=output_file, width=7.2, height=11.5, units="in")



cell_type_name = "Adipocytes"
output_file <- paste0(visualization_dir, tissue_10_model_stem, cell_type_name, "_loadings_scatter_colored_by_sex.pdf")
scatters <- make_cell_type_loadings_scatters_no_gtex_color(loadings, covariates$Adipocytes, cell_type_name, covariates$sex, "sex")
ggsave(scatters, file=output_file, width=7.2, height=10.5, units="in")

cell_type_name = "Neurons"
output_file <- paste0(visualization_dir, tissue_10_model_stem, cell_type_name, "_loadings_scatter.pdf")
scatters <- make_cell_type_loadings_scatters(loadings, covariates$Neurons, cell_type_name, tissue_10_names, tissue_colors)
ggsave(scatters, file=output_file, width=7.2, height=10.5, units="in")

cell_type_name = "Myocytes"
output_file <- paste0(visualization_dir, tissue_10_model_stem, cell_type_name, "_loadings_scatter.pdf")
scatters <- make_cell_type_loadings_scatters(loadings, covariates$Myocytes, cell_type_name, tissue_10_names, tissue_colors)
ggsave(scatters, file=output_file, width=7.2, height=10.5, units="in")


cell_type_name = "Keratinocytes"
output_file <- paste0(visualization_dir, tissue_10_model_stem, cell_type_name, "_loadings_scatter.pdf")
scatters <- make_cell_type_loadings_scatters(loadings, covariates$Keratinocytes, cell_type_name, tissue_10_names, tissue_colors)
ggsave(scatters, file=output_file, width=7.2, height=10.5, units="in")


cell_type_name = "Epithelial_cells"
output_file <- paste0(visualization_dir, tissue_10_model_stem, cell_type_name, "_loadings_scatter.pdf")
scatters <- make_cell_type_loadings_scatters(loadings, covariates$Epithelial, cell_type_name, tissue_10_names, tissue_colors)
ggsave(scatters, file=output_file, width=7.2, height=10.5, units="in")

cell_type_name = "Adipocytes"
output_file <- paste0(visualization_dir, tissue_10_model_stem, cell_type_name, "_loadings_scatter.pdf")
scatters <- make_cell_type_loadings_scatters(loadings, covariates$Adipocytes, cell_type_name, tissue_10_names, tissue_colors)
ggsave(scatters, file=output_file, width=7.2, height=10.5, units="in")

cell_type_name = "Neutrophils"
output_file <- paste0(visualization_dir, tissue_10_model_stem, cell_type_name, "_loadings_scatter.pdf")
scatters <- make_cell_type_loadings_scatters(loadings, covariates$Neutrophils, cell_type_name, tissue_10_names, tissue_colors)
ggsave(scatters, file=output_file, width=7.2, height=10.5, units="in")


cell_type_name = "Age"
output_file <- paste0(visualization_dir, tissue_10_model_stem, cell_type_name, "_loadings_scatter.pdf")
scatters <- make_cell_type_loadings_scatters(loadings, covariates$age, cell_type_name, tissue_10_names, tissue_colors)
ggsave(scatters, file=output_file, width=7.2, height=10.5, units="in")


cell_type_name = "ischemic_time"
output_file <- paste0(visualization_dir, tissue_10_model_stem, cell_type_name, "_loadings_scatter.pdf")
scatters <- make_cell_type_loadings_scatters(loadings, covariates$ischemic_time, cell_type_name, tissue_10_names, tissue_colors)
ggsave(scatters, file=output_file, width=7.2, height=10.5, units="in")

cell_type_name = "bmi"
output_file <- paste0(visualization_dir, tissue_10_model_stem, cell_type_name, "_loadings_scatter.pdf")
scatters <- make_cell_type_loadings_scatters(loadings, covariates$bmi, cell_type_name, tissue_10_names, tissue_colors)
ggsave(scatters, file=output_file, width=7.2, height=10.5, units="in")

cell_type_name = "weight"
output_file <- paste0(visualization_dir, tissue_10_model_stem, cell_type_name, "_loadings_scatter.pdf")
scatters <- make_cell_type_loadings_scatters(loadings, covariates$weight, cell_type_name, tissue_10_names, tissue_colors)
ggsave(scatters, file=output_file, width=7.2, height=10.5, units="in")


cell_type_name = "height"
output_file <- paste0(visualization_dir, tissue_10_model_stem, cell_type_name, "_loadings_scatter.pdf")
scatters <- make_cell_type_loadings_scatters(loadings, covariates$height, cell_type_name, tissue_10_names, tissue_colors)
ggsave(scatters, file=output_file, width=7.2, height=10.5, units="in")


cell_type_name = "genotype_pc1"
output_file <- paste0(visualization_dir, tissue_10_model_stem, cell_type_name, "_loadings_scatter.pdf")
scatters <- make_cell_type_loadings_scatters(loadings, pcs$genotype_PC0, cell_type_name, tissue_10_names, tissue_colors)
ggsave(scatters, file=output_file, width=7.2, height=10.5, units="in")


cell_type_name = "genotype_pc2"
output_file <- paste0(visualization_dir, tissue_10_model_stem, cell_type_name, "_loadings_scatter.pdf")
scatters <- make_cell_type_loadings_scatters(loadings, pcs$genotype_PC1, cell_type_name, tissue_10_names, tissue_colors)
ggsave(scatters, file=output_file, width=7.2, height=10.5, units="in")



















}
















