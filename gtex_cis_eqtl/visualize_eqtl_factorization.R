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

make_loading_boxplot_plot_by_race <- function(sample_covariate_file, loading_file) {
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
		#race_vec <- c(race_vec, as.character(as.numeric(covariates$race==3)))
		race_vec <- c(race_vec, covariates$race)
		factor_vec <- c(factor_vec, rep(as.character(factor_number), length(covariates$race)))
	}


	df <- data.frame(loading=loading_vec, race=factor(race_vec), latent_factor=factor(factor_vec, levels=as.character(1:num_factors)))

	boxplot <- ggplot(df, aes(x=latent_factor, y=loading, fill=race)) + geom_boxplot(outlier.size = .1) +
				gtex_v8_figure_theme() + 
	        	labs(x="Latent factor", y = "Sample loading", fill="Known race") +
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

make_loading_boxplot_plot_by_sex <- function(sample_covariate_file, loading_file) {
	#tissues <- read.table(tissue_file, header=FALSE)
	loadings <- read.table(loading_file, header=FALSE)
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

make_loading_boxplot_plot_by_cohort <- function(sample_covariate_file, loading_file) {
	#tissues <- read.table(tissue_file, header=FALSE)
	loadings <- read.table(loading_file, header=FALSE)
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

make_loading_boxplot_plot_by_tissue <- function(tissues,tissue_colors, loading_file) {
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
	        	labs(x="Latent factor", y = "Sample loading", fill="") +
	        	theme(legend.position="bottom") +
	           	guides(fill=guide_legend(nrow=5,byrow=TRUE, override.aes = list(size=.1))) + 
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
	           labs(x=loading_name, y=cell_type_name, title=paste0(tissue_name, " / correlation pvalue: ", corry$p.value))

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

explore_relationship_between_surveyed_covariates_and_eqtl_factors <- function(surveyed_covariate_file, loading_file) {
	loadings <- read.table(loading_file, header=FALSE)
	covs <- read.table(surveyed_covariate_file, header=TRUE)
	covs <- covs[,2:(dim(covs)[2])]
	factor_num <- 11
	num_covs = dim(covs)[2]
	for (cov_num in 1:num_covs) {
		cov_vec <- covs[, cov_num]
		loading_vec <- loadings[,factor_num]

		valid_samples <- cov_vec == 0.0 | cov_vec == 1
		lin_model <- lm(loading_vec[valid_samples] ~ cov_vec[valid_samples])
		print(paste0(cov_num, "    ", summary(lin_model)$adj.r.squared))
	}
	#print(summary(covs))
}


explore_relationship_between_technical_covariates_and_eqtl_factors <- function(surveyed_covariate_file, loading_file, indi_names) {
	loadings <- read.table(loading_file, header=FALSE)
	covs <- read.table(surveyed_covariate_file, header=TRUE, sep="\t")
	covs <- covs[,2:(dim(covs)[2])]
	print(colnames(covs))
	factor_num <- 7
	num_covs = dim(covs)[2]
	loading_vec <- loadings[, factor_num]
	indices <- loading_vec > .1
	cov_vec <- covs[, 3]
	indi_names <- factor(indi_names)

	lin_model <- lm(loading_vec ~ indi_names)
	print(summary(lin_model))

	for (cov_num in 1:num_covs) {
		cov_vec <- covs[, cov_num]
		loading_vec <- loadings[,factor_num]

		lin_model <- lm(loading_vec ~ cov_vec)
		print(paste0(cov_num, "    ", summary(lin_model)$adj.r.squared))
	}
	#print(summary(covs))
}



processed_data_dir <- args[1]
eqtl_results_dir <- args[2]
visualization_dir <- args[3]
tissue_colors_file <- args[4]


# Read in tissue colors and names
tissue_colors = read.table(tissue_colors_file, header = T, stringsAsFactors = F, sep = "\t")
# slight mislabeling
for (tiss_num in 1:length(tissue_colors$tissue_id)) {
	if (tissue_colors$tissue_id[tiss_num] == "Brain_Spinal_cord_cervical_c1") {
		tissue_colors$tissue_id[tiss_num] = "Brain_Spinal_cord_cervical_c.1"
	}
	if (tissue_colors$tissue_id[tiss_num] == "Cells_EBVtransformed_lymphocytes") {
		tissue_colors$tissue_id[tiss_num] = "Cells_EBV.transformed_lymphocytes"
	}
}


############################
# Load in files
############################
tissue_10_file <- paste0(processed_data_dir, "tissues_subset_10_sample_names.txt")
tissue_10_sample_covariate_file <- paste0(processed_data_dir, "tissues_subset_10_sample_covariates.txt")
tissue_10_surveyed_covariate_file <- paste0(processed_data_dir, "tissues_subset_10_sample_surveyed_covariates.txt")
tissue_10_technical_covariate_file <- paste0(processed_data_dir, "tissues_subset_10_sample_technical_covariates.txt")

tissue_10_names <- get_tissue_names(tissue_10_file)
tissue_10_indi_names <- get_indi_names(tissue_10_file)


############################
# Model Specification
############################
tissue_10_model_stem <- paste0("tissues_subset_10_lf_interaction_egenes_eqtl_factorization_vi_results_k_init_20_lambda_v_1_seed_6_RE_init2_temper_")
tissue_10_loading_file <- paste0(eqtl_results_dir, tissue_10_model_stem, "U_S.txt")



############################
# Start making plots!!
############################

#explore_relationship_between_surveyed_covariates_and_eqtl_factors(tissue_10_surveyed_covariate_file, tissue_10_loading_file)
explore_relationship_between_technical_covariates_and_eqtl_factors(tissue_10_technical_covariate_file, tissue_10_loading_file, tissue_10_indi_names)

######################
# Make box plot for each Race, showing loading distributions
output_file <- paste0(visualization_dir, tissue_10_model_stem, "race_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_race(tissue_10_sample_covariate_file, tissue_10_loading_file)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")


######################
# Make box plot for each Cohort, showing loading distributions
output_file <- paste0(visualization_dir, tissue_10_model_stem, "cohort_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_cohort(tissue_10_sample_covariate_file, tissue_10_loading_file)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")

######################
# Make box plot for each Sex, showing loading distributions
output_file <- paste0(visualization_dir, tissue_10_model_stem, "sex_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_sex(tissue_10_sample_covariate_file, tissue_10_loading_file)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")

######################
# Make box plot for each tissue, showing loading distributions
output_file <- paste0(visualization_dir, tissue_10_model_stem, "tissue_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_tissue(tissue_10_names, tissue_colors, tissue_10_loading_file)
print(output_file)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")


##################
# Make scatter plot showing correlation between loadings and cell type enrichments
covariates <- read.table(tissue_10_sample_covariate_file, header=TRUE, sep="\t")
loadings <- read.table(tissue_10_loading_file, header=FALSE)


if (FALSE) {
# Epithelial cells
output_file <- paste0(visualization_dir, tissue_10_model_stem, "epiethelial_cells_loadings_scatter.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatters(loadings, covariates$Epithelial_cells, "Epithelial", tissue_10_names, tissue_colors)
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=11, units="in")


# Epithelial cells
output_file <- paste0(visualization_dir, tissue_10_model_stem, "keratinocytes_loadings_scatter.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatters(loadings, covariates$Keratinocytes, "Keratinocytes", tissue_10_names, tissue_colors)
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=11, units="in")


# Myocyte cells
output_file <- paste0(visualization_dir, tissue_10_model_stem, "myocyte_loadings_scatter.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatters(loadings, covariates$Myocytes, "Myocytes", tissue_10_names, tissue_colors)
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=11, units="in")
}

if (FALSE) {

temp_cov_file <- "/work-zfs/abattle4/bstrober/qtl_factorization/gtex_cis_eqtl/processed_data/tissues_subset_10_residual_expression_covariates.txt"
temp_cov <- read.table(temp_cov_file, header=TRUE, sep="\t")


valid_tissues <- c("Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg")
loading_num <- 1
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_1_genotype_PC1_scatter_colored_by_skin_tissues_and_race.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_colored_by_specified_tissues_and_categorical_variables(loadings[, loading_num], temp_cov$genotype_PC0, "Genotype PC1", "Skin Tissue", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, factor(covariates$race==2, levels=c(TRUE, FALSE)), "African Ancestry")
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")


loading_num <- 1
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_1_genotype_PC1_scatter_colored_by_tissue_and_race.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_colored_by_tissues_and_categorical_variables(loadings[, loading_num], temp_cov$genotype_PC0, "Genotype PC1", paste0("Loading ", loading_num), tissue_10_names, tissue_colors, valid_tissues, factor(covariates$race==2, levels=c(TRUE, FALSE)), "African Ancestry")
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")


valid_tissues <- c("Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg")
loading_num <- 1
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_1_genotype_PC1_scatter_for_skin_tissues_colored_by_race.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue_colored_by_categorical_variable(loadings[, loading_num], temp_cov$genotype_PC0, "Genotype PC1", "Skin Tissue", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, factor(covariates$race), "Ancestry")
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")
}

if (FALSE) {
##################
# Make scatter plot showing correlation between 2 loadings and color samples by cell type enrichments for samples from subset of tissues
valid_tissues <- c("Thyroid")
loading_num1 <- 3
loading_num2 <- 6
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_3_loading_6_scatter_colored_by_Epithelial_for_thyroid_samples.pdf")
ct_loading_scatter <- make_loadings_loadings_scatter_colored_by_cell_type_for_samples_from_specified_tissue(loadings[, loading_num1], loadings[, loading_num2],  covariates$Epithelial, "Epithelial", "Thyroid tissue", paste0("Loading ", loading_num1), paste0("Loading ", loading_num2), tissue_10_names, valid_tissues)
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")

valid_tissues <- c("Colon_Sigmoid", "Small_Intestine_Terminal_Ileum", "Stomach")
loading_num1 <- 4
loading_num2 <- 7
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_4_loading_7_scatter_colored_by_Epithelial_for_digestive_samples.pdf")
ct_loading_scatter <- make_loadings_loadings_scatter_colored_by_cell_type_for_samples_from_specified_tissue(loadings[, loading_num1], loadings[, loading_num2],  covariates$Epithelial, "Epithelial", "Digestive tissues", paste0("Loading ", loading_num1), paste0("Loading ", loading_num2), tissue_10_names, valid_tissues)
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")


##################
# Make scatter plot showing correlation between loadings and cell type enrichments for samples from subset of tissues
valid_tissues <- c("Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg")
loading_num <- 1
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_1_Epithelial_scatter_for_skin_samples.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue(loadings[, loading_num], covariates$Epithelial, "Epithelial", "Skin tissues", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, tissue_colors)
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")

valid_tissues <- c("Colon_Sigmoid", "Small_Intestine_Terminal_Ileum", "Stomach")
loading_num <- 4
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_4_Epithelial_scatter_for_digestive_samples.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue(loadings[, loading_num], covariates$Epithelial, "Epithelial", "Digestive tissues", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, tissue_colors)
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")

valid_tissues <- c("Colon_Sigmoid", "Small_Intestine_Terminal_Ileum", "Stomach")
loading_num <- 7
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_7_Epithelial_scatter_for_digestive_samples.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue(loadings[, loading_num], covariates$Epithelial, "Epithelial", "Digestive tissues", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, tissue_colors)
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")


valid_tissues <- c("Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg")
loading_num <- 1
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_1_keratinocytes_scatter_for_skin_samples.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue(loadings[, loading_num], covariates$Keratinocytes, "Keratinocytes", "Skin tissues", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, tissue_colors)
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")

valid_tissues <- c("Esophagus_Mucosa")
loading_num <- 5
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_5_Epithelial_scatter_for_esophagus_mucosa.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue(loadings[, loading_num], covariates$Epithelial, "Epithelial", "Esophagus_Mucosa", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, tissue_colors)
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")

valid_tissues <- c("Esophagus_Mucosa")
loading_num <- 5
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_5_keratinocytes_scatter_for_esophagus_mucosa.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue(loadings[, loading_num], covariates$Keratinocytes, "Keratinocytes", "Esophagus_Mucosa", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, tissue_colors)
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")

valid_tissues <- c("Pituitary")
loading_num <- 9
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_9_Epithelial_scatter_for_pituitary.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue(loadings[, loading_num], covariates$Epithelial, "Epithelial", "Pituitary", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, tissue_colors)
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")

valid_tissues <- c("Thyroid")
loading_num <- 3
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_3_Epithelial_scatter_for_thyroid.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue(loadings[, loading_num], covariates$Epithelial, "Epithelial", "Thyroid", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, tissue_colors)
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")

valid_tissues <- c("Thyroid")
loading_num <- 6
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_6_Epithelial_scatter_for_thyroid.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue(loadings[, loading_num], covariates$Epithelial, "Epithelial", "Thyroid", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, tissue_colors)
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")


valid_tissues <- c("Thyroid")
loading_num <- 3
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_3_Ischemic_time_scatter_for_thyroid.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue(loadings[, loading_num], covariates$ischemic_time, "ischemic time", "Thyroid", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, tissue_colors)
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")

valid_tissues <- c("Thyroid")
loading_num <- 6
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_6_Ischemic_time_scatter_for_thyroid.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue(loadings[, loading_num], covariates$ischemic_time, "ischemic time", "Thyroid", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, tissue_colors)
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")




##################
# Make scatter plot showing correlation between loadings and cell type enrichments for samples from subset of tissues colored by cohort
valid_tissues <- c("Pituitary")
loading_num <- 9
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_9_Epithelial_scatter_for_pituitary_colored_by_cohort.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue_colored_by_categorical_variable(loadings[, loading_num], covariates$Epithelial, "Epithelial", "Pituitary", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, covariates$cohort, "")
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")

valid_tissues <- c("Esophagus_Mucosa")
loading_num <- 5
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_5_keratinocytes_scatter_for_esophagus_mucosa_colored_by_cohort.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue_colored_by_categorical_variable(loadings[, loading_num], covariates$Keratinocytes, "Keratinocytes", "Esophagus_Mucosa", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, covariates$cohort, "")
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")

valid_tissues <- c("Esophagus_Mucosa")
loading_num <- 5
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_5_Epithelial_scatter_for_esophagus_mucosa_colored_by_cohort.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue_colored_by_categorical_variable(loadings[, loading_num], covariates$Epithelial, "Epithelial", "Esophagus_Mucosa", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, covariates$cohort, "")
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")


valid_tissues <- c("Colon_Sigmoid", "Small_Intestine_Terminal_Ileum", "Stomach")
loading_num <- 7
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_7_Epithelial_scatter_for_digestive_samples_colored_by_cohort.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue_colored_by_categorical_variable(loadings[, loading_num], covariates$Epithelial, "Epithelial", "Digestive tissues", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, covariates$cohort, "")
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")

valid_tissues <- c("Colon_Sigmoid", "Small_Intestine_Terminal_Ileum", "Stomach")
loading_num <- 4
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_4_Epithelial_scatter_for_digestive_samples_colored_by_cohort.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue_colored_by_categorical_variable(loadings[, loading_num], covariates$Epithelial, "Epithelial", "Digestive tissues", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, covariates$cohort, "")
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")


valid_tissues <- c("Thyroid")
loading_num <- 3
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_3_Epithelial_scatter_for_thyroid_colored_by_cohort.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue_colored_by_categorical_variable(loadings[, loading_num], covariates$Epithelial, "Epithelial", "Thyroid", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, covariates$cohort, "")
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")

valid_tissues <- c("Thyroid")
loading_num <- 6
output_file <- paste0(visualization_dir, tissue_10_model_stem, "loading_6_Epithelial_scatter_for_thyroid_colored_by_cohort.pdf")
ct_loading_scatter <- make_cell_type_loadings_scatter_for_samples_from_specified_tissue_colored_by_categorical_variable(loadings[, loading_num], covariates$Epithelial, "Epithelial", "Thyroid", paste0("Loading ", loading_num), tissue_10_names, valid_tissues, covariates$cohort, "")
ggsave(ct_loading_scatter, file=output_file, width=7.2, height=5, units="in")
}

if (FALSE) {
#####################
# Run Umap on loadings. Plot Umap loadings in scatter plot color by observed Xcell cell type enrichments
loadings <- read.table(tissue_10_loading_file, header=FALSE)
umap_loadings = umap(loadings)$layout
covariates <- read.table(tissue_10_sample_covariate_file, header=TRUE, sep="\t")

# Adipocytes
output_file <- paste0(visualization_dir, tissue_10_model_stem, "umap_loadings_scatter_colored_by_adipocyte_enrichments.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_vector(umap_loadings, covariates$Adipocytes, "Adipocyte enrichment")
ggsave(umap_scatter, file=output_file, width=7.2, height=5.5, units="in")

# Epithelial_cells
output_file <- paste0(visualization_dir, tissue_10_model_stem, "umap_loadings_scatter_colored_by_epithelial_cell_enrichments.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_vector(umap_loadings, covariates$Epithelial_cells, "Epithelial cells enrichment")
ggsave(umap_scatter, file=output_file, width=7.2, height=5.5, units="in")


# Hepatocytes
output_file <- paste0(visualization_dir, tissue_10_model_stem, "umap_loadings_scatter_colored_by_hepatocytes_cell_enrichments.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_vector(umap_loadings, covariates$Hepatocytes, "Hepatocytes enrichment")
ggsave(umap_scatter, file=output_file, width=7.2, height=5.5, units="in")


# Keratinocytes
output_file <- paste0(visualization_dir, tissue_10_model_stem, "umap_loadings_scatter_colored_by_keratinocytes_cell_enrichments.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_vector(umap_loadings, covariates$Keratinocytes, "Keratinocytes enrichment")
ggsave(umap_scatter, file=output_file, width=7.2, height=5.5, units="in")


# Myocytes
output_file <- paste0(visualization_dir, tissue_10_model_stem, "umap_loadings_scatter_colored_by_myocytes_cell_enrichments.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_vector(umap_loadings, covariates$Myocytes, "Myocytes enrichment")
ggsave(umap_scatter, file=output_file, width=7.2, height=5.5, units="in")

# Neurons
output_file <- paste0(visualization_dir, tissue_10_model_stem, "umap_loadings_scatter_colored_by_neurons_cell_enrichments.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_vector(umap_loadings, covariates$Neurons, "Neurons enrichment")
ggsave(umap_scatter, file=output_file, width=7.2, height=5.5, units="in")


# Neutrophils
output_file <- paste0(visualization_dir, tissue_10_model_stem, "umap_loadings_scatter_colored_by_neutrophils_cell_enrichments.pdf")
umap_scatter <- make_umap_loading_scatter_plot_colored_by_real_valued_vector(umap_loadings, covariates$Neutrophils, "Neutrophils enrichment")
ggsave(umap_scatter, file=output_file, width=7.2, height=5.5, units="in")


#####################
# Run Umap on loadings. Plot Umap loadings in scatter plot color by observed tissue type
output_file <- paste0(visualization_dir, tissue_10_model_stem, "umap_loading_scatter.pdf")
umap_scatter <- make_umap_loading_scatter_plot(tissue_10_names, tissue_colors, tissue_10_sample_covariate_file, tissue_10_loading_file, umap_loadings)
ggsave(umap_scatter, file=output_file, width=7.2*1.5, height=5.5*1.5, units="in")

}











































###############
## OLD
###############







if (FALSE) {
######################
# Make box plot for each Race, showing loading distributions
# 4-tissue (No SVI)
output_file <- paste0(visualization_dir, tissue_4_no_svi_model_stem, "race_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_race(tissue_4_sample_covariate_file, tissue_4_no_svi_loading_file)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")
# 4 tissue (SVI)
output_file <- paste0(visualization_dir, tissue_4_svi_model_stem, "race_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_race(tissue_4_sample_covariate_file, tissue_4_svi_loading_file)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")
# 20 tissue (SVI)
output_file <- paste0(visualization_dir, tissue_20_model_stem, "race_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_race(tissue_20_sample_covariate_file, tissue_20_loading_file)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")


######################
# Make box plot for each tissue, showing loading distributions
# 4-tissue (No SVI)
output_file <- paste0(visualization_dir, tissue_4_no_svi_model_stem, "tissue_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_tissue(tissue_4_names, tissue_colors, tissue_4_no_svi_loading_file)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")
# 4-tissue (SVI)
output_file <- paste0(visualization_dir, tissue_4_svi_model_stem, "tissue_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_tissue(tissue_4_names, tissue_colors, tissue_4_svi_loading_file)
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")
# 20-tissue (SVI)
output_file <- paste0(visualization_dir, tissue_20_model_stem, "tissue_colored_loading_boxplot.pdf")
boxplot <- make_loading_boxplot_plot_by_tissue(tissue_20_names, tissue_colors, tissue_20_loading_file)
ggsave(boxplot, file=output_file, width=14, height=5.5, units="in")



#####################
# Run Umap on loadings. Plot Umap loadings in scatter plot color by observed tissue type
# 4-tissue (No SVI)
output_file <- paste0(visualization_dir, tissue_4_no_svi_model_stem, "umap_loading_scatter.pdf")
umap_scatter <- make_umap_loading_scatter_plot(tissue_4_names, tissue_colors, tissue_4_sample_covariate_file, tissue_4_no_svi_loading_file)
ggsave(umap_scatter, file=output_file, width=7.2*1.5, height=5.5*1.5, units="in")
# 4-tissue (SVI)
output_file <- paste0(visualization_dir, tissue_4_svi_model_stem, "umap_loading_scatter.pdf")
umap_scatter <- make_umap_loading_scatter_plot(tissue_4_names, tissue_colors, tissue_4_sample_covariate_file, tissue_4_svi_loading_file)
ggsave(umap_scatter, file=output_file, width=7.2*1.5, height=5.5*1.5, units="in")
# 20-tissue (SVI)
output_file <- paste0(visualization_dir, tissue_20_model_stem, "umap_loading_scatter.pdf")
umap_scatter <- make_umap_loading_scatter_plot(tissue_20_names, tissue_colors, tissue_20_sample_covariate_file, tissue_20_loading_file)
ggsave(umap_scatter, file=output_file, width=7.2*1.5, height=5.5*1.5, units="in")


######################
# Make scatter plot where each sample is a point, x and y axis are factor loadings, and points are colored by their tissue type
# 4-tissue (No SVI)
output_file <- paste0(visualization_dir, tissue_4_no_svi_model_stem, "loading_scatter.pdf")
scatter <- make_loading_scatter_plot(tissue_4_names, tissue_4_sample_covariate_file, tissue_colors, tissue_4_no_svi_loading_file)
ggsave(scatter, file=output_file, width=7.2, height=5.5, units="in")
# 4-tissue (SVI)
output_file <- paste0(visualization_dir, tissue_4_svi_model_stem, "loading_scatter.pdf")
scatter <- make_loading_scatter_plot(tissue_4_names, tissue_4_sample_covariate_file, tissue_colors, tissue_4_svi_loading_file)
ggsave(scatter, file=output_file, width=7.2, height=5.5, units="in")




######################
# Make heatmap comparing (showing correlation of) loading matrices of two models trained on the same data
output_file <- paste0(visualization_dir, "gtex_4_tissue_svi_vs_no_svi_loading_correlation_heatmap.pdf")
svi_vs_no_svi_heatmap <- make_loading_correlation_heatmap(tissue_4_no_svi_loading_file, tissue_4_svi_loading_file, "Loadings (standard VI)", "Loadings (Stochastic VI)")
ggsave(svi_vs_no_svi_heatmap, file=output_file, width=7.2, height=5.5, units="in")



######################
# Make box-plot of factor scores of "ancestry-factor" stratefied by whether there is missingness in only 1 of ancestry studies
# 4-tissue (No SVI)
output_file <- paste0(visualization_dir, tissue_4_no_svi_model_stem, "ancestry_factor_boxplot_colored_by_delta_ancestry_eqtl_missingness.pdf")
ancestry_boxplot <- make_ancestry_missingness_boxplot(tissue_4_no_svi_factor_file, unique(tissue_4_names), tissue_4_ancestry_specific_eqtl_file_root)
ggsave(ancestry_boxplot, file=output_file, width=7.2, height=5.0, units="in")

######################
# Make Scatter-plot showing median-delta (european ancestry eqtl beta vs african ancestry eqtl beta) across tissues against by factor scores of "ancestry-factor"
# 4-tissue (No SVI)
output_file <- paste0(visualization_dir, tissue_4_no_svi_model_stem, "median_delta_ancestry_specific_effect_size_scatter_vs_ancestry_factor_score.pdf")
ancestry_scatter <- make_median_ancestry_scatter_plot(tissue_4_no_svi_factor_file, unique(tissue_4_names), tissue_4_ancestry_specific_eqtl_file_root)
ggsave(ancestry_scatter, file=output_file, width=7.2, height=5.5, units="in")

######################
# Make Scatter-plot showing delta (european ancestry eqtl beta vs african ancestry eqtl beta) against by factor scores of "ancestry-factor"
# 4-tissue (No SVI)
output_file <- paste0(visualization_dir, tissue_4_no_svi_model_stem, "delta_ancestry_specific_effect_size_scatter_vs_ancestry_factor_score.pdf")
ancestry_scatter <- make_ancestry_scatter_plot(tissue_4_no_svi_factor_file, unique(tissue_4_names), tissue_4_ancestry_specific_eqtl_file_root)
ggsave(ancestry_scatter, file=output_file, width=7.2, height=8.5, units="in")
}





