args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(cowplot)
library(umap)
library(ggplot2)
library(RColorBrewer)
library(sigmoid)
library(lme4)
library(plyr)
options(bitmapType = 'cairo', device = 'pdf')



gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


extract_num_coloc_hits_standard_vs_surge_df <- function(study_names, coloc_dir) {
	gwas_study_name_vec <- c()
	eqtl_study_name_vec <- c()
	num_hits_vec <- c()
	pph4_thresh_vec <- c()

	pph4_threshs <- c(.8, .9, .95, .99)

	for (pph4_iter in 1:length(pph4_threshs)) {
		pph4_thresh = pph4_threshs[pph4_iter]
		for (study_iter in 1:length(study_names)) {
			gwas_study_name <- study_names[study_iter]
			surge_genes <- c()
			for (latent_factor_num in 1:10) {
				eqtl_study_name <- paste0("surge_latent_factor_", latent_factor_num, "_interaction")
				results_file <- paste0(coloc_dir, eqtl_study_name, "_", gwas_study_name, "_coloc_test_results.txt")
				coloc_results_df <- read.table(results_file, header=TRUE)
				coloc_genes <- as.character(coloc_results_df[coloc_results_df$pph4 > pph4_thresh,]$gene_name)
				surge_genes <- c(surge_genes, coloc_genes)
			}
			num_hits <- length(unique(surge_genes))
			gwas_study_name_vec <- c(gwas_study_name_vec, gwas_study_name)
			eqtl_study_name_vec <- c(eqtl_study_name_vec, "SURGE_interaction_eqtl")
			num_hits_vec <- c(num_hits_vec, num_hits)
			pph4_thresh_vec <- c(pph4_thresh_vec, pph4_thresh)
		}
	}

	eqtl_study_name <- "standard_eqtl"
	informal_eqtl_study_name <- "standard_eqtl"
	for (pph4_iter in 1:length(pph4_threshs)) {
		pph4_thresh = pph4_threshs[pph4_iter]
		for (study_iter in 1:length(study_names)) {
			gwas_study_name <- study_names[study_iter]
			results_file <- paste0(coloc_dir, eqtl_study_name, "_", gwas_study_name, "_coloc_test_results.txt")
			coloc_results_df <- read.table(results_file, header=TRUE)
			num_hits <- sum(coloc_results_df$pph4 > pph4_thresh)

			gwas_study_name_vec <- c(gwas_study_name_vec, gwas_study_name)
			eqtl_study_name_vec <- c(eqtl_study_name_vec, informal_eqtl_study_name)
			num_hits_vec <- c(num_hits_vec, num_hits)
			pph4_thresh_vec <- c(pph4_thresh_vec, pph4_thresh)
		}
	}

	df <- data.frame(gwas_study=factor(gwas_study_name_vec), eqtl_study=factor(eqtl_study_name_vec), num_colocalizations=num_hits_vec, pph4_thresh=pph4_thresh_vec)
	return(df)
}



extract_num_coloc_hits_df <- function(study_names, coloc_dir) {
	gwas_study_name_vec <- c()
	eqtl_study_name_vec <- c()
	num_hits_vec <- c()
	pph4_thresh_vec <- c()

	pph4_threshs <- c(.8, .9, .95, .99)

	for (latent_factor_num in 1:10) {
		eqtl_study_name <- paste0("surge_latent_factor_", latent_factor_num, "_interaction")
		informal_eqtl_study_name <- paste0("SURGE_interaction_", latent_factor_num)
		for (pph4_iter in 1:length(pph4_threshs)) {
			pph4_thresh = pph4_threshs[pph4_iter]
			for (study_iter in 1:length(study_names)) {
				gwas_study_name <- study_names[study_iter]
				results_file <- paste0(coloc_dir, eqtl_study_name, "_", gwas_study_name, "_coloc_test_results.txt")
				coloc_results_df <- read.table(results_file, header=TRUE)
				num_hits <- sum(coloc_results_df$pph4 > pph4_thresh)

				gwas_study_name_vec <- c(gwas_study_name_vec, gwas_study_name)
				eqtl_study_name_vec <- c(eqtl_study_name_vec, informal_eqtl_study_name)
				num_hits_vec <- c(num_hits_vec, num_hits)
				pph4_thresh_vec <- c(pph4_thresh_vec, pph4_thresh)
			}
		}
	}

	eqtl_study_name <- "standard_eqtl"
	informal_eqtl_study_name <- "standard_eqtl"
	for (pph4_iter in 1:length(pph4_threshs)) {
		pph4_thresh = pph4_threshs[pph4_iter]
		for (study_iter in 1:length(study_names)) {
			gwas_study_name <- study_names[study_iter]
			results_file <- paste0(coloc_dir, eqtl_study_name, "_", gwas_study_name, "_coloc_test_results.txt")
			coloc_results_df <- read.table(results_file, header=TRUE)
			num_hits <- sum(coloc_results_df$pph4 > pph4_thresh)

			gwas_study_name_vec <- c(gwas_study_name_vec, gwas_study_name)
			eqtl_study_name_vec <- c(eqtl_study_name_vec, informal_eqtl_study_name)
			num_hits_vec <- c(num_hits_vec, num_hits)
			pph4_thresh_vec <- c(pph4_thresh_vec, pph4_thresh)
		}
	}

	df <- data.frame(gwas_study=factor(gwas_study_name_vec), eqtl_study=factor(eqtl_study_name_vec), num_colocalizations=num_hits_vec, pph4_thresh=pph4_thresh_vec)
	return(df)
}

make_number_of_colocalizations_bar_plot <- function(df, pph4_threshold) {
	df$gwas_study <- factor(df$gwas_study, levels=c("ukbb_bmi", "ukbb_eczema", "sle", "ukbb_blood_eosinophil_count", "ukbb_blood_high_light_scatter_reticulotye_count", "ukbb_blood_lymphocyte_count", "ukbb_blood_mean_corpuscular_hemoglobin", "ukbb_blood_monocyte_count", "ukbb_blood_platelet_count", "ukbb_blood_platelet_vol", "ukbb_blood_red_count", "ukbb_blood_white_count"))
	df$gwas_study = revalue(df$gwas_study, c("ukbb_bmi"="bmi", "ukbb_eczema"="eczema", "ukbb_blood_eosinophil_count"="eosinophil count", "ukbb_blood_high_light_scatter_reticulotye_count"="reticulocyte count", "ukbb_blood_lymphocyte_count"="lymphocyte count", "ukbb_blood_mean_corpuscular_hemoglobin"="corpuscular hemoglobin", "ukbb_blood_monocyte_count"="monocyte count", "ukbb_blood_platelet_count"="platelet count", "ukbb_blood_red_count"="blood red count", "ukbb_blood_white_count"="blood white count"))
	p <- ggplot(df, aes(fill=eqtl_study, y=num_colocalizations, x=gwas_study)) + 
    	geom_bar(position="dodge", stat="identity") +
    	gtex_v8_figure_theme() + theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
    	theme(legend.position="top") +
    	labs(y="Number of colocalizations", x="GWAS study", fill="", title=paste0("pph4 threshold: ", pph4_threshold))
    return(p)
}


processed_gwas_studies_file <- args[1]
coloc_dir <- args[2]
visualize_coloc_dir <- args[3]


study_df <- read.table(processed_gwas_studies_file, header=FALSE)
study_names <- as.character(study_df$V1)


num_coloc_hits_standard_vs_surge_df <- extract_num_coloc_hits_standard_vs_surge_df(study_names, coloc_dir)

num_coloc_hits_df <- extract_num_coloc_hits_df(study_names, coloc_dir)


pph4_threshold <- .9
output_file <- paste0(visualize_coloc_dir, "number_of_colocalizations_pph4_thresh_", pph4_threshold, "_bar_plot.pdf")
bar_plot_9 <- make_number_of_colocalizations_bar_plot(num_coloc_hits_df[num_coloc_hits_df$pph4_thresh==pph4_threshold,], pph4_threshold)
ggsave(bar_plot_9, file=output_file, width=7.2, height=6.0, units="in")

pph4_threshold <- .95
output_file <- paste0(visualize_coloc_dir, "number_of_colocalizations_pph4_thresh_", pph4_threshold, "_bar_plot.pdf")
bar_plot_95 <- make_number_of_colocalizations_bar_plot(num_coloc_hits_df[num_coloc_hits_df$pph4_thresh==pph4_threshold,], pph4_threshold)
ggsave(bar_plot_95, file=output_file, width=7.2, height=6.0, units="in")

pph4_threshold <- .99
output_file <- paste0(visualize_coloc_dir, "number_of_colocalizations_pph4_thresh_", pph4_threshold, "_bar_plot.pdf")
bar_plot_99 <- make_number_of_colocalizations_bar_plot(num_coloc_hits_df[num_coloc_hits_df$pph4_thresh==pph4_threshold,], pph4_threshold)
ggsave(bar_plot_99, file=output_file, width=7.2, height=6.0, units="in")






pph4_threshold <- .9
output_file <- paste0(visualize_coloc_dir, "number_of_colocalizations_standard_v_surge_pph4_thresh_", pph4_threshold, "_bar_plot.pdf")
bar_plot_9 <- make_number_of_colocalizations_bar_plot(num_coloc_hits_standard_vs_surge_df[num_coloc_hits_standard_vs_surge_df$pph4_thresh==pph4_threshold,], pph4_threshold)
ggsave(bar_plot_9, file=output_file, width=10.2, height=6.0, units="in")

pph4_threshold <- .95
output_file <- paste0(visualize_coloc_dir, "number_of_colocalizations_standard_v_surge_pph4_thresh_", pph4_threshold, "_bar_plot.pdf")
bar_plot_95 <- make_number_of_colocalizations_bar_plot(num_coloc_hits_standard_vs_surge_df[num_coloc_hits_standard_vs_surge_df$pph4_thresh==pph4_threshold,], pph4_threshold)
ggsave(bar_plot_95, file=output_file, width=10.2, height=6.0, units="in")

pph4_threshold <- .99
output_file <- paste0(visualize_coloc_dir, "number_of_colocalizations_standard_v_surge_pph4_thresh_", pph4_threshold, "_bar_plot.pdf")
bar_plot_99 <- make_number_of_colocalizations_bar_plot(num_coloc_hits_standard_vs_surge_df[num_coloc_hits_standard_vs_surge_df$pph4_thresh==pph4_threshold,], pph4_threshold)
ggsave(bar_plot_99, file=output_file, width=10.2, height=6.0, units="in")




