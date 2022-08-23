args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(cowplot)
library(umap)
library(ggplot2)
library(RColorBrewer)
library(lme4)
options(bitmapType = 'cairo', device = 'pdf')



gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}



extract_alphabetical_ordered_gene_level_pvalues_across_factors <- function(stem, num_factors) {
	factor_num <- 1

	pvalues <- list()

	for (factor_num in 1:num_factors) {

		file_name <- paste0(stem, factor_num, "_alphabetical_ordered_gene_level_stats.txt")

		temp_data <- read.table(file_name, header=TRUE)

		aaa = -log10(temp_data$pvalue*temp_data$num_snps_in_gene + 1e-300)
		print(length(aaa))

		pvalues[[factor_num]] <- aaa

	}

	return(pvalues)

}


make_real_perm_correlation_matrix_at_gene_level <- function(real_list, perm_list, num_factors) {
	corr_mat <- matrix(0, num_factors, num_factors)

	for (factor1 in 1:10) {
		for (factor2 in 1:10) {
			if (factor1 != factor2) {
				corry = cor(real_list[[factor1]], perm_list[[factor2]])
				corr_mat[factor1, factor2] = corry
			}
		}
	}

    melted_mat <- melt(corr_mat)
    colnames(melted_mat) <- c("Real", "Permutation", "Correlation")




  	melted_mat$Real = factor(melted_mat$Real, levels=as.character(1:num_factors))
  	melted_mat$Permutation = factor(melted_mat$Permutation, levels=as.character(1:num_factors))

    
    #levels(melted_mat$PC) = paste0("PC", 1:(length(levels(melted_mat$PC))))
    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=Real, y=Permutation)) + geom_tile(aes(fill=Correlation)) + scale_fill_gradient2(midpoint=.0, guide="colorbar") + 
    gtex_v8_figure_theme() +
    theme(text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
    # Save File
    return(heatmap)
}


efdr_calculation <- function(real_pvalues, perm_pvalues) {
	sorted_real = sort(real_pvalues)
	M1 = length(real_pvalues)
	M2 = length(perm_pvalues)
	eFDRs <- c()
	prev_min = 10000
	for (index in 1:length(sorted_real)) {
		real_value = sorted_real[index]
		frac_real = sum(sorted_real >= real_value)/M1
		frac_perm = sum(perm_pvalues >= real_value)/M2

		curr_efdr = frac_perm/frac_real
		if (curr_efdr > prev_min) {
			curr_efdr = prev_min
		}
		prev_min = curr_efdr
		eFDRs <- c(eFDRs, curr_efdr)
	}
	return(eFDRs)
}

make_real_vs_perm_gene_qq_plot <- function(real_pvalues, perm_pvalues, latent_factor) {
	efdrs = efdr_calculation(real_pvalues, perm_pvalues)
	num_egenes_05 <- sum(efdrs < .05)
	num_egenes_2 <- sum(efdrs < .2)

	df <- data.frame(real=sort(real_pvalues), perm=sort(perm_pvalues))

	 plotter <- ggplot(df) + 
             geom_point(aes(x=perm, y=real), size=.1) +
             geom_abline()+
             gtex_v8_figure_theme() + 
             labs(x="Permutation", y = "Real", title=paste0("Context ", latent_factor, "\n", num_egenes_05, " eGenes (eFDR < .05)\n", num_egenes_2, " eGenes (eFDR < .2)")) + 
             theme(legend.text = element_text(size=8), legend.title = element_text(size=8)) 
  return(plotter)

}


gene_qq_plot_colored_by_real_or_perm <- function(real_gene_level_pvalues, perm_gene_level_pvalues, num_factors) {
	num_genes <- length(real_gene_level_pvalues[[1]])
	null_pvalues = runif(num_genes)
	real <- c()
	null <- c()
	context <- c()
	version <- c()
	ordered_contexts <-c()

	for (lf in 1:num_factors) {
		sorted_null = sort(-log10(null_pvalues + 1e-300))
		sorted_real = sort(real_gene_level_pvalues[[lf]])

		real <- c(real, sorted_real)
		null <- c(null, sorted_null)
		context <- c(context, rep(paste0("context_", lf), num_genes))
		version <- c(version, rep("real", num_genes))



		sorted_null = sort(-log10(null_pvalues + 1e-300))
		sorted_real = sort(perm_gene_level_pvalues[[lf]])

		real <- c(real, sorted_real)
		null <- c(null, sorted_null)
		context <- c(context, rep(paste0("context_", lf), num_genes))
		version <- c(version, rep("perm", num_genes))


		ordered_contexts <- c(ordered_contexts, paste0("context_", lf))
	}
	df <- data.frame(real=real, null=null, version=version)
	 plotter <- ggplot(df) + 
             geom_point(aes(x=null, y=real, color=version), size=1) +
             geom_abline()+
             gtex_v8_figure_theme() + 
             labs(x="Permutation", y = "Real", color="") + 
             theme(legend.text = element_text(size=8), legend.title = element_text(size=8)) 
    return(plotter)

}


gene_qq_plot_colored_by_factors <- function(gene_level_pvalues, num_factors, pvalue_type) {


	num_genes <- length(gene_level_pvalues[[1]])
	null_pvalues = runif(num_genes)

	real <- c()
	null <- c()
	context <- c()

	ordered_contexts <-c()
	for (lf in 1:num_factors) {
		sorted_null = sort(-log10(null_pvalues + 1e-300))
		sorted_real = sort(gene_level_pvalues[[lf]])

		real <- c(real, sorted_real)
		null <- c(null, sorted_null)
		context <- c(context, rep(paste0("context_", lf), num_genes))
		ordered_contexts <- c(ordered_contexts, paste0("context_", lf))
	}
	df <- data.frame(real=real, null=null, context=factor(context,levels=ordered_contexts))
	 plotter <- ggplot(df) + 
             geom_point(aes(x=null, y=real, color=context), size=1) +
             geom_abline()+
             gtex_v8_figure_theme() + 
             labs(x="Permutation", y = "Real", title=pvalue_type) + 
             theme(legend.text = element_text(size=8), legend.title = element_text(size=8)) 
    return(plotter)
}


output_stem <- args[1]

num_factors = 1

real_gene_level_pvalues_alphabetical_list <- extract_alphabetical_ordered_gene_level_pvalues_across_factors(paste0(output_stem, "False_interaction_eqtl_results_latent_factor_"), num_factors)
perm_gene_level_pvalues_alphabetical_list <- extract_alphabetical_ordered_gene_level_pvalues_across_factors(paste0(output_stem, "interaction_only_interaction_eqtl_results_latent_factor_"), num_factors)

gene_qq_plots <- list()
union_perm_gene_level_pvalues <- c()
for (latent_factor in 1:1) {
	union_perm_gene_level_pvalues <- c(union_perm_gene_level_pvalues, perm_gene_level_pvalues_alphabetical_list[[latent_factor]])
}



for (latent_factor in 1:1) {
	efdrs = efdr_calculation(real_gene_level_pvalues_alphabetical_list[[latent_factor]], union_perm_gene_level_pvalues)
	num_egenes_05 <- sum(efdrs < .05)
	num_egenes_2 <- sum(efdrs < .2)
}




output_file <- paste0(output_stem, "real_and_perm_vs_null_gene_qq_plot.pdf")
#real_and_perm_vs_null_gene_qq_plot <- gene_qq_plot_colored_by_real_or_perm(real_gene_level_pvalues_alphabetical_list, perm_gene_level_pvalues_alphabetical_list, num_factors)
#ggsave(real_and_perm_vs_null_gene_qq_plot, file=output_file, width=7.2, height=6, units="in")


output_file <- paste0(output_stem, "real_vs_null_gene_qq_plot.pdf")
#real_null_gene_qq_plot <- gene_qq_plot_colored_by_factors(real_gene_level_pvalues_alphabetical_list, num_factors, "Real")
#ggsave(real_null_gene_qq_plot, file=output_file, width=7.2, height=6, units="in")

output_file <- paste0(output_stem, "perm_vs_null_gene_qq_plot.pdf")
#perm_null_gene_qq_plot <- gene_qq_plot_colored_by_factors(perm_gene_level_pvalues_alphabetical_list, num_factors, "Permuted")
#ggsave(perm_null_gene_qq_plot, file=output_file, width=7.2, height=6, units="in")



output_file <- paste0(output_stem, "real_perm_gene_level_correlation_heatmap.pdf")
#real_perm_corr_heatmap <- make_real_perm_correlation_matrix_at_gene_level(real_gene_level_pvalues_alphabetical_list, perm_gene_level_pvalues_alphabetical_list, num_factors)
#ggsave(real_perm_corr_heatmap, file=output_file, width=7.2, height=6, units="in")

gene_qq_plots <- list()
for (latent_factor in 1:1) {
	output_file <- paste0(output_stem, "gene_qq_plot_real_vs_perm_", latent_factor, ".pdf")
	gene_qq_plots[[latent_factor]] <- make_real_vs_perm_gene_qq_plot(real_gene_level_pvalues_alphabetical_list[[latent_factor]], perm_gene_level_pvalues_alphabetical_list[[latent_factor]], latent_factor)
	ggsave(gene_qq_plots[[latent_factor]], file=output_file, width=7.2, height=6, units="in")
}
output_file <- paste0(output_stem, "gene_qq_plot_real_vs_perm_all_factors.pdf")
#merged_qq = plot_grid(plotlist=gene_qq_plots, ncol=3)
#ggsave(merged_qq, file=output_file, width=7.2, height=9.3, units="in")




