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

	pph4_threshs <- c(.5, .6, .7, .8, .9, .95, .99)

	for (pph4_iter in 1:length(pph4_threshs)) {
		pph4_thresh = pph4_threshs[pph4_iter]
		for (study_iter in 1:length(study_names)) {
			gwas_study_name <- study_names[study_iter]
			surge_genes <- c()
			for (latent_factor_num in 1:6) {
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

	for (pph4_iter in 1:length(pph4_threshs)) {
		pph4_thresh = pph4_threshs[pph4_iter]
		for (study_iter in 1:length(study_names)) {
			gwas_study_name <- study_names[study_iter]
			surge_genes <- c()
			for (latent_factor_num in 1:6) {
				eqtl_study_name <- paste0("expression_pc_", latent_factor_num, "_interaction")
				results_file <- paste0(coloc_dir, eqtl_study_name, "_", gwas_study_name, "_coloc_test_results.txt")
				coloc_results_df <- read.table(results_file, header=TRUE)
				coloc_genes <- as.character(coloc_results_df[coloc_results_df$pph4 > pph4_thresh,]$gene_name)
				surge_genes <- c(surge_genes, coloc_genes)
			}

			num_hits <- length(unique(surge_genes))
			#num_hits <- length((surge_genes))
			gwas_study_name_vec <- c(gwas_study_name_vec, gwas_study_name)
			eqtl_study_name_vec <- c(eqtl_study_name_vec, "expression_pc_interaction_eqtl")
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

extract_num_coloc_hits_standard_vs_surge_df_tirtiary <- function(study_names, coloc_dir) {
	gwas_study_name_vec <- c()
	eqtl_study_name_vec <- c()
	num_hits_vec <- c()
	pph4_thresh_vec <- c()

	pph4_threshs <- c(.5, .6, .7, .8, .9, .95, .99)

	latent_factors <- c(1, 2, 4)

	for (pph4_iter in 1:length(pph4_threshs)) {
		pph4_thresh = pph4_threshs[pph4_iter]
		for (study_iter in 1:length(study_names)) {
			gwas_study_name <- study_names[study_iter]
			surge_genes <- c()
			for (latent_factor_iter in 1:length(latent_factors)) {
				latent_factor_num = latent_factors[latent_factor_iter]
				eqtl_study_name <- paste0("surge_latent_factor_", latent_factor_num, "_interaction")
				results_file <- paste0(coloc_dir, eqtl_study_name, "_", gwas_study_name, "_coloc_test_results.txt")
				coloc_results_df <- read.table(results_file, header=TRUE)
				coloc_genes <- as.character(coloc_results_df[coloc_results_df$pph4 > pph4_thresh,]$gene_name)
				surge_genes <- c(surge_genes, coloc_genes)
			}

			num_hits <- length(unique(surge_genes))
			#num_hits <- length((surge_genes))
			gwas_study_name_vec <- c(gwas_study_name_vec, gwas_study_name)
			eqtl_study_name_vec <- c(eqtl_study_name_vec, "SURGE_interaction_eqtl")
			num_hits_vec <- c(num_hits_vec, num_hits)
			pph4_thresh_vec <- c(pph4_thresh_vec, pph4_thresh)
		}
	}

	latent_factors <- c(1, 2, 3)
	for (pph4_iter in 1:length(pph4_threshs)) {
		pph4_thresh = pph4_threshs[pph4_iter]
		for (study_iter in 1:length(study_names)) {
			gwas_study_name <- study_names[study_iter]
			surge_genes <- c()
			for (latent_factor_iter in 1:length(latent_factors)) {
				latent_factor_num = latent_factors[latent_factor_iter]
				eqtl_study_name <- paste0("expression_pc_", latent_factor_num, "_interaction")
				results_file <- paste0(coloc_dir, eqtl_study_name, "_", gwas_study_name, "_coloc_test_results.txt")
				coloc_results_df <- read.table(results_file, header=TRUE)
				coloc_genes <- as.character(coloc_results_df[coloc_results_df$pph4 > pph4_thresh,]$gene_name)
				surge_genes <- c(surge_genes, coloc_genes)
			}

			num_hits <- length(unique(surge_genes))
			#num_hits <- length((surge_genes))
			gwas_study_name_vec <- c(gwas_study_name_vec, gwas_study_name)
			eqtl_study_name_vec <- c(eqtl_study_name_vec, "expression_pc_interaction_eqtl")
			num_hits_vec <- c(num_hits_vec, num_hits)
			pph4_thresh_vec <- c(pph4_thresh_vec, pph4_thresh)
		}
	}

	if (FALSE) {
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
	}

	df <- data.frame(gwas_study=factor(gwas_study_name_vec), eqtl_study=factor(eqtl_study_name_vec), num_colocalizations=num_hits_vec, pph4_thresh=pph4_thresh_vec)
	

	return(df)
}

extract_num_coloc_hits_standard_vs_surge_df_secondary <- function(study_names, coloc_dir) {
	gwas_study_name_vec <- c()
	eqtl_study_name_vec <- c()
	num_hits_vec <- c()
	pph4_thresh_vec <- c()

	pph4_threshs <- c(.5,.6, .7, .8, .9, .95, .99)

	latent_factors <- c(3, 5, 6)
	#latent_factors <- c(1, 2, 4)

	for (pph4_iter in 1:length(pph4_threshs)) {
		pph4_thresh = pph4_threshs[pph4_iter]
		for (study_iter in 1:length(study_names)) {
			gwas_study_name <- study_names[study_iter]
			surge_genes <- c()
			for (latent_factor_iter in 1:length(latent_factors)) {
				latent_factor_num = latent_factors[latent_factor_iter]
				eqtl_study_name <- paste0("surge_latent_factor_", latent_factor_num, "_interaction")
				results_file <- paste0(coloc_dir, eqtl_study_name, "_", gwas_study_name, "_coloc_test_results.txt")
				coloc_results_df <- read.table(results_file, header=TRUE)
				coloc_genes <- as.character(coloc_results_df[coloc_results_df$pph4 > pph4_thresh,]$gene_name)
				surge_genes <- c(surge_genes, coloc_genes)
			}

			num_hits <- length(unique(surge_genes))
			#num_hits <- length((surge_genes))
			gwas_study_name_vec <- c(gwas_study_name_vec, gwas_study_name)
			eqtl_study_name_vec <- c(eqtl_study_name_vec, "SURGE_interaction_eqtl")
			num_hits_vec <- c(num_hits_vec, num_hits)
			pph4_thresh_vec <- c(pph4_thresh_vec, pph4_thresh)
		}
	}

	#latent_factors <- c(1, 2, 3)
	latent_factors <- c(4, 5, 6)
	for (pph4_iter in 1:length(pph4_threshs)) {
		pph4_thresh = pph4_threshs[pph4_iter]
		for (study_iter in 1:length(study_names)) {
			gwas_study_name <- study_names[study_iter]
			surge_genes <- c()
			for (latent_factor_iter in 1:length(latent_factors)) {
				latent_factor_num = latent_factors[latent_factor_iter]
				eqtl_study_name <- paste0("expression_pc_", latent_factor_num, "_interaction")
				results_file <- paste0(coloc_dir, eqtl_study_name, "_", gwas_study_name, "_coloc_test_results.txt")
				coloc_results_df <- read.table(results_file, header=TRUE)
				coloc_genes <- as.character(coloc_results_df[coloc_results_df$pph4 > pph4_thresh,]$gene_name)
				surge_genes <- c(surge_genes, coloc_genes)
			}

			num_hits <- length(unique(surge_genes))
			#num_hits <- length((surge_genes))
			gwas_study_name_vec <- c(gwas_study_name_vec, gwas_study_name)
			eqtl_study_name_vec <- c(eqtl_study_name_vec, "expression_pc_interaction_eqtl")
			num_hits_vec <- c(num_hits_vec, num_hits)
			pph4_thresh_vec <- c(pph4_thresh_vec, pph4_thresh)
		}
	}

	if (FALSE) {
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

	for (latent_factor_num in 1:6) {
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

	for (latent_factor_num in 1:6) {
		eqtl_study_name <- paste0("expression_pc_", latent_factor_num, "_interaction")
		informal_eqtl_study_name <- paste0("expression_pc_", latent_factor_num)
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
	print('#####################')
	print(pph4_threshold)
	#print(df$num_colocalizations[as.character(df$eqtl_study) == "expression_pc_interaction_eqtl"])
	#print(df$num_colocalizations[as.character(df$eqtl_study) == "SURGE_interaction_eqtl"])

	#print(wilcox.test(df$num_colocalizations[as.character(df$eqtl_study) == "expression_pc_interaction_eqtl"], df$num_colocalizations[as.character(df$eqtl_study) == "SURGE_interaction_eqtl"], paired=TRUE,alternative = "greater"))
	#print(wilcox.test(df$num_colocalizations[as.character(df$eqtl_study) == "SURGE_interaction_eqtl"], df$num_colocalizations[as.character(df$eqtl_study) == "expression_pc_interaction_eqtl"],  paired=TRUE,alternative = "greater"))
	print(t.test(df$num_colocalizations[as.character(df$eqtl_study) == "SURGE_interaction_eqtl"], df$num_colocalizations[as.character(df$eqtl_study) == "expression_pc_interaction_eqtl"], paired=TRUE))
	print(wilcox.test(df$num_colocalizations[as.character(df$eqtl_study) == "SURGE_interaction_eqtl"], df$num_colocalizations[as.character(df$eqtl_study) == "expression_pc_interaction_eqtl"], paired=TRUE))
	print(t.test(df$num_colocalizations[as.character(df$eqtl_study) == "SURGE_interaction_eqtl"], df$num_colocalizations[as.character(df$eqtl_study) == "standard_eqtl"], paired=TRUE))
	print(wilcox.test(df$num_colocalizations[as.character(df$eqtl_study) == "SURGE_interaction_eqtl"], df$num_colocalizations[as.character(df$eqtl_study) == "standard_eqtl"], paired=TRUE))



	#df$gwas_study <- factor(df$gwas_study, levels=c("ukbb_bmi", "ukbb_eczema", "sle", "ukbb_blood_eosinophil_count", "ukbb_blood_high_light_scatter_reticulotye_count", "ukbb_blood_lymphocyte_count", "ukbb_blood_mean_corpuscular_hemoglobin", "ukbb_blood_monocyte_count", "ukbb_blood_platelet_count", "ukbb_blood_platelet_vol", "ukbb_blood_red_count", "ukbb_blood_white_count"))
	#df$gwas_study = revalue(df$gwas_study, c("ukbb_bmi"="bmi", "ukbb_eczema"="eczema", "ukbb_blood_eosinophil_count"="eosinophil count", "ukbb_blood_high_light_scatter_reticulotye_count"="reticulocyte count", "ukbb_blood_lymphocyte_count"="lymphocyte count", "ukbb_blood_mean_corpuscular_hemoglobin"="corpuscular hemoglobin", "ukbb_blood_monocyte_count"="monocyte count", "ukbb_blood_platelet_count"="platelet count", "ukbb_blood_red_count"="blood red count", "ukbb_blood_white_count"="blood white count"))
	df$gwas_study = revalue(df$gwas_study, c("ukbb_blood_mean_corpuscular_hemoglobin"="corpuscular hemoglobin", "ukbb_blood_monocyte_count"="monocyte count","ukbb_blood_platelet_vol"="platelet volume", "ukbb_blood_high_light_scatter_reticulotye_count"="reticulocyte count", "ukbb_eczema"="eczema", "ukbb_cholesterol"="cholesterol", "ukbb_height"="height", "ukbb_bone_mineral_density"="BMD", "ukbb_hair_color"="hair color", "ukbb_FEV1FVCz"="FEV1FVC", "ukbb_balding"="balding", "ukbb_fvc"="FVC", "ukbb_diastolic_bp"="diastolic BP", "ukbb_menarche_age"="menarche age", "ukbb_number_children"="no. children"))
	df$gwas_study = factor(df$gwas_study, levels=c("corpuscular hemoglobin", "monocyte count","platelet volume", "reticulocyte count", "eczema", "cholesterol", "height", "BMD", "hair color", "FEV1FVC", "balding", "FVC", "diastolic BP", "menarche age", "no. children"))

	df$eqtl_study = factor(df$eqtl_study, levels=c("standard_eqtl", "expression_pc_interaction_eqtl", "SURGE_interaction_eqtl"))


	p <- ggplot(df, aes(fill=eqtl_study, y=num_colocalizations, x=gwas_study)) + 
    	geom_bar(position="dodge", stat="identity") +
    	gtex_v8_figure_theme() + theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
    	theme(legend.position="top") +
    	labs(y="Number of colocalizations", x="GWAS study", fill="", title=paste0("pph4 threshold: ", pph4_threshold))
    return(p)
}

make_number_of_colocalizations_bar_plot2 <- function(df, pph4_threshold) {

	#print(df$num_colocalizations[as.character(df$eqtl_study) == "standard_eqtl"])
	#print(df$num_colocalizations[as.character(df$eqtl_study) == "SURGE_interaction_eqtl"])

	#print(wilcox.test(df$num_colocalizations[as.character(df$eqtl_study) == "expression_pc_interaction_eqtl"], df$num_colocalizations[as.character(df$eqtl_study) == "SURGE_interaction_eqtl"], paired=TRUE,alternative = "greater"))
	#print(wilcox.test(df$num_colocalizations[as.character(df$eqtl_study) == "SURGE_interaction_eqtl"], df$num_colocalizations[as.character(df$eqtl_study) == "expression_pc_interaction_eqtl"],  paired=TRUE,alternative = "greater"))
	#print(t.test(df$num_colocalizations[as.character(df$eqtl_study) == "SURGE_interaction_eqtl"], df$num_colocalizations[as.character(df$eqtl_study) == "standard_eqtl"], paired=TRUE))

	df = df[startsWith(as.character(df$eqtl_study), "expression") == FALSE,]

	#df$gwas_study <- factor(df$gwas_study, levels=c("ukbb_bmi", "ukbb_eczema", "sle", "ukbb_blood_eosinophil_count", "ukbb_blood_high_light_scatter_reticulotye_count", "ukbb_blood_lymphocyte_count", "ukbb_blood_mean_corpuscular_hemoglobin", "ukbb_blood_monocyte_count", "ukbb_blood_platelet_count", "ukbb_blood_platelet_vol", "ukbb_blood_red_count", "ukbb_blood_white_count"))
	#df$gwas_study = revalue(df$gwas_study, c("ukbb_bmi"="bmi", "ukbb_eczema"="eczema", "ukbb_blood_eosinophil_count"="eosinophil count", "ukbb_blood_high_light_scatter_reticulotye_count"="reticulocyte count", "ukbb_blood_lymphocyte_count"="lymphocyte count", "ukbb_blood_mean_corpuscular_hemoglobin"="corpuscular hemoglobin", "ukbb_blood_monocyte_count"="monocyte count", "ukbb_blood_platelet_count"="platelet count", "ukbb_blood_red_count"="blood red count", "ukbb_blood_white_count"="blood white count"))
	df$gwas_study = revalue(df$gwas_study, c("ukbb_blood_mean_corpuscular_hemoglobin"="corpuscular hemoglobin", "ukbb_blood_monocyte_count"="monocyte count","ukbb_blood_platelet_vol"="platelet volume", "ukbb_blood_high_light_scatter_reticulotye_count"="reticulocyte count", "ukbb_eczema"="eczema", "ukbb_cholesterol"="cholesterol", "ukbb_height"="height", "ukbb_bone_mineral_density"="BMD", "ukbb_hair_color"="hair color", "ukbb_FEV1FVCz"="FEV1FVC", "ukbb_balding"="balding", "ukbb_fvc"="FVC", "ukbb_diastolic_bp"="diastolic BP", "ukbb_menarche_age"="menarche age", "ukbb_number_children"="no. children"))
	df$gwas_study = factor(df$gwas_study, levels=c("corpuscular hemoglobin", "monocyte count","platelet volume", "reticulocyte count", "eczema", "cholesterol", "height", "BMD", "hair color", "FEV1FVC", "balding", "FVC", "diastolic BP", "menarche age", "no. children"))

	df$eqtl_study = revalue(df$eqtl_study, c("standard_eqtl"="standard_eQTL", "SURGE_interaction_1"="SURGE_1_eQTL","SURGE_interaction_2"="SURGE_2_eQTL", "SURGE_interaction_3"="SURGE_3_eQTL", "SURGE_interaction_4"="SURGE_4_eQTL", "SURGE_interaction_5"="SURGE_5_eQTL", "SURGE_interaction_6"="SURGE_6_eQTL"))

	#df$eqtl_study = factor(df$eqtl_study, levels=c("standard_eqtl", "SURGE_interaction_eqtl"))


	p <- ggplot(df, aes(fill=eqtl_study, y=num_colocalizations, x=gwas_study)) + 
    	geom_bar(position="dodge", stat="identity") +
    	gtex_v8_figure_theme() + theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
    	theme(legend.position="top") +
    	labs(y="Number of colocalizations", x="GWAS study", fill="", title=paste0("pph4 threshold: ", pph4_threshold))
    return(p)
}

make_number_of_coloc_plot_for_single_component <- function(df, pph4_threshold, component_number) {
	df2 = df[as.character(df$eqtl_study) == paste0("expression_pc_", component_number) | as.character(df$eqtl_study) == paste0("SURGE_", component_number, "_eQTL") ,]

	# Gosh I suck at R
	if (component_number == 1) {
		df2$eqtl_study = revalue(df2$eqtl_study, c("expression_pc_1"="expression pc", "SURGE_1_eQTL"="SURGE"))
	}
	if (component_number == 2) {
		df2$eqtl_study = revalue(df2$eqtl_study, c("expression_pc_2"="expression pc", "SURGE_2_eQTL"="SURGE"))
	}
	if (component_number == 3) {
		df2$eqtl_study = revalue(df2$eqtl_study, c("expression_pc_3"="expression pc", "SURGE_3_eQTL"="SURGE"))
	}
	if (component_number == 4) {
		df2$eqtl_study = revalue(df2$eqtl_study, c("expression_pc_4"="expression pc", "SURGE_4_eQTL"="SURGE"))
	}
	if (component_number == 5) {
		df2$eqtl_study = revalue(df2$eqtl_study, c("expression_pc_5"="expression pc", "SURGE_5_eQTL"="SURGE"))
	}
	if (component_number == 6) {
		df2$eqtl_study = revalue(df2$eqtl_study, c("expression_pc_6"="expression pc", "SURGE_6_eQTL"="SURGE"))
	}


	df2$eqtl_study = factor(df2$eqtl_study, levels=c("expression pc", "SURGE"))

	print(df2[as.character(df2$gwas_study)=="monocyte count",])

	p <- ggplot(df2, aes(fill=eqtl_study, y=num_colocalizations, x=gwas_study)) + 
    	geom_bar(position="dodge", stat="identity") +
    	gtex_v8_figure_theme() + theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
    	theme(legend.position="top") +
    	scale_fill_manual(values = c("grey", "seagreen3")) +
    	labs(y="No.\n colocalizations", x="GWAS study", fill="", title=paste0("Component/Context ", component_number))
    return(p)
}


make_number_of_colocalizations_bar_stratefied_by_component <- function(df, pph4_threshold) {
	df = df[startsWith(as.character(df$eqtl_study), "standard") == FALSE,]

	df$gwas_study = revalue(df$gwas_study, c("ukbb_blood_mean_corpuscular_hemoglobin"="corpuscular hemoglobin", "ukbb_blood_monocyte_count"="monocyte count","ukbb_blood_platelet_vol"="platelet volume", "ukbb_blood_high_light_scatter_reticulotye_count"="reticulocyte count", "ukbb_eczema"="eczema", "ukbb_cholesterol"="cholesterol", "ukbb_height"="height", "ukbb_bone_mineral_density"="BMD", "ukbb_hair_color"="hair color", "ukbb_FEV1FVCz"="FEV1FVC", "ukbb_balding"="balding", "ukbb_fvc"="FVC", "ukbb_diastolic_bp"="diastolic BP", "ukbb_menarche_age"="menarche age", "ukbb_number_children"="no. children"))
	df$gwas_study = factor(df$gwas_study, levels=c("corpuscular hemoglobin", "monocyte count","platelet volume", "reticulocyte count", "eczema", "cholesterol", "height", "BMD", "hair color", "FEV1FVC", "balding", "FVC", "diastolic BP", "menarche age", "no. children"))

	df$eqtl_study = revalue(df$eqtl_study, c("standard_eqtl"="standard_eQTL", "SURGE_interaction_1"="SURGE_1_eQTL","SURGE_interaction_2"="SURGE_2_eQTL", "SURGE_interaction_3"="SURGE_3_eQTL", "SURGE_interaction_4"="SURGE_4_eQTL", "SURGE_interaction_5"="SURGE_5_eQTL", "SURGE_interaction_6"="SURGE_6_eQTL"))

	component_number = 1
	p1 <- make_number_of_coloc_plot_for_single_component(df, pph4_threshold, component_number)+ theme(axis.text.x=element_blank()) + labs(x="")
	legender <- get_legend(p1)
	p1 <- p1 + theme(legend.position="none")

	component_number = 2
	p2 <- make_number_of_coloc_plot_for_single_component(df, pph4_threshold, component_number)+ theme(legend.position="none") + theme(axis.text.x=element_blank())+ labs(x="")

	component_number = 3
	p3 <- make_number_of_coloc_plot_for_single_component(df, pph4_threshold, component_number)+ theme(legend.position="none")+ theme(axis.text.x=element_blank())+ labs(x="")

	component_number = 4
	p4 <- make_number_of_coloc_plot_for_single_component(df, pph4_threshold, component_number)+ theme(legend.position="none")+ theme(axis.text.x=element_blank())+ labs(x="")

	component_number = 5
	p5 <- make_number_of_coloc_plot_for_single_component(df, pph4_threshold, component_number)+ theme(legend.position="none")+ theme(axis.text.x=element_blank())+ labs(x="")

	component_number = 6
	p6 <- make_number_of_coloc_plot_for_single_component(df, pph4_threshold, component_number)+ theme(legend.position="none")

	joint_plot <- plot_grid(legender, p1, p2, p3, p4, p5, p6, ncol=1, rel_heights=c(.3, 1, 1, 1, 1, 1, 2.3))


	return(joint_plot)
}

make_surge_expression_pc_coloc_comparison_scatter <- function(df, pph4_threshold) {
	df <- df[(df$max_expression_pc_pph4 > pph4_threshold) | (df$max_surge_pph4 > pph4_threshold),]
	p <- ggplot(df, aes(x=max_expression_pc_pph4, y=max_surge_pph4, color=blood_trait)) + geom_point() +
	gtex_v8_figure_theme() +
	labs(x="Coloc PPH4 from expression-PC eQTL", y="Coloc PPH4 from expression-PC eQT")
	return(p)
}

make_surge_expression_pc_coloc_comparison_density <- function(df, pph4_threshold) {
	df <- df[(df$max_expression_pc_pph4 > pph4_threshold) | (df$max_surge_pph4 > pph4_threshold),]

	pph4s <- c(df$max_expression_pc_pph4, df$max_surge_pph4)
	eqtl_class <- c(rep("expression_pc", length(df$max_expression_pc_pph4)), rep("SURGE", length(df$max_surge_pph4)))

	df2 <- data.frame(pph4=pph4s, eqtl_class=factor(eqtl_class))

	p <- ggplot(df2, aes(pph4, colour = eqtl_class)) +
  		geom_density() +
  		gtex_v8_figure_theme()

	return(p)
}

processed_gwas_studies_file <- args[1]
coloc_dir <- args[2]
visualize_coloc_dir <- args[3]

study_df <- read.table(processed_gwas_studies_file, header=FALSE)
study_names <- as.character(study_df$V1)
study_names <- study_names[1:15]
print(study_names)


# Scatter plot comparring surge and expression pcs
pph4_threshold = .9
surge_expr_pc_cmp_df <- read.table(paste0(coloc_dir, "surge_vs_expression_pc_coloc_hits.txt"), header=TRUE, sep="\t")
output_file <- paste0(visualize_coloc_dir, "SURGE_expr_pc_coloc_", pph4_threshold, "_comparison_scatter.pdf")
scatter <- make_surge_expression_pc_coloc_comparison_scatter(surge_expr_pc_cmp_df, pph4_threshold)
ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")


# Scatter plot comparring surge and expression pcs
pph4_threshold = .1
surge_expr_pc_cmp_df <- read.table(paste0(coloc_dir, "surge_vs_expression_pc_coloc_hits.txt"), header=TRUE, sep="\t")
output_file <- paste0(visualize_coloc_dir, "SURGE_expr_pc_coloc_", pph4_threshold, "_comparison_density.pdf")
density <- make_surge_expression_pc_coloc_comparison_density(surge_expr_pc_cmp_df, pph4_threshold)
ggsave(density, file=output_file, width=7.2, height=3.0, units="in")


if (FALSE) {
num_coloc_hits_standard_vs_surge_df <- extract_num_coloc_hits_standard_vs_surge_df(study_names, coloc_dir)
save(num_coloc_hits_standard_vs_surge_df,file="num_coloc_hits_standard_vs_surge_df.Rda")
num_coloc_hits_unshared_standard_vs_surge_df <- extract_num_coloc_hits_standard_vs_surge_df_secondary(study_names, coloc_dir)
save(num_coloc_hits_unshared_standard_vs_surge_df,file="num_coloc_hits_unshared_standard_vs_surge_df.Rda")
num_coloc_hits_shared_standard_vs_surge_df <- extract_num_coloc_hits_standard_vs_surge_df_tirtiary(study_names, coloc_dir)
save(num_coloc_hits_shared_standard_vs_surge_df,file="num_coloc_hits_shared_standard_vs_surge_df.Rda")
}
load("num_coloc_hits_standard_vs_surge_df.Rda")
load("num_coloc_hits_unshared_standard_vs_surge_df.Rda")
load("num_coloc_hits_shared_standard_vs_surge_df.Rda")



pph4_threshold <- .5
output_file <- paste0(visualize_coloc_dir, "number_of_colocalizations_standard_v_surge_pph4_thresh_", pph4_threshold, "_bar_plot.pdf")
bar_plot_9 <- make_number_of_colocalizations_bar_plot(num_coloc_hits_standard_vs_surge_df[num_coloc_hits_standard_vs_surge_df$pph4_thresh==pph4_threshold,], pph4_threshold)
ggsave(bar_plot_9, file=output_file, width=10.2, height=6.0, units="in")

pph4_threshold <- .7
output_file <- paste0(visualize_coloc_dir, "number_of_colocalizations_standard_v_surge_pph4_thresh_", pph4_threshold, "_bar_plot.pdf")
bar_plot_9 <- make_number_of_colocalizations_bar_plot(num_coloc_hits_standard_vs_surge_df[num_coloc_hits_standard_vs_surge_df$pph4_thresh==pph4_threshold,], pph4_threshold)
ggsave(bar_plot_9, file=output_file, width=10.2, height=6.0, units="in")


pph4_threshold <- .8
output_file <- paste0(visualize_coloc_dir, "number_of_colocalizations_standard_v_surge_pph4_thresh_", pph4_threshold, "_bar_plot.pdf")
bar_plot_9 <- make_number_of_colocalizations_bar_plot(num_coloc_hits_standard_vs_surge_df[num_coloc_hits_standard_vs_surge_df$pph4_thresh==pph4_threshold,], pph4_threshold)
ggsave(bar_plot_9, file=output_file, width=10.2, height=6.0, units="in")

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








if (FALSE) {
num_coloc_hits_per_component_df <- extract_num_coloc_hits_df(study_names, coloc_dir)
save(num_coloc_hits_per_component_df,file="num_coloc_hits_per_component_df.Rda")
}
load("num_coloc_hits_per_component_df.Rda")

pph4_threshold <- .9
output_file <- paste0(visualize_coloc_dir, "number_of_colocalizations_pph4_thresh_", pph4_threshold, "_bar_plot.pdf")
bar_plot_9 <- make_number_of_colocalizations_bar_plot2(num_coloc_hits_per_component_df[num_coloc_hits_per_component_df$pph4_thresh==pph4_threshold,], pph4_threshold)
ggsave(bar_plot_9, file=output_file, width=7.2, height=6.0, units="in")

pph4_threshold <- .95
output_file <- paste0(visualize_coloc_dir, "number_of_colocalizations_pph4_thresh_", pph4_threshold, "_bar_plot.pdf")
bar_plot_95 <- make_number_of_colocalizations_bar_plot2(num_coloc_hits_per_component_df[num_coloc_hits_per_component_df$pph4_thresh==pph4_threshold,], pph4_threshold)
ggsave(bar_plot_95, file=output_file, width=7.2, height=6.0, units="in")

pph4_threshold <- .99
output_file <- paste0(visualize_coloc_dir, "number_of_colocalizations_pph4_thresh_", pph4_threshold, "_bar_plot.pdf")
bar_plot_99 <- make_number_of_colocalizations_bar_plot2(num_coloc_hits_per_component_df[num_coloc_hits_per_component_df$pph4_thresh==pph4_threshold,], pph4_threshold)
ggsave(bar_plot_99, file=output_file, width=7.2, height=6.0, units="in")



pph4_threshold <- .95
joint_bar_plot <- make_number_of_colocalizations_bar_stratefied_by_component(num_coloc_hits_per_component_df[num_coloc_hits_per_component_df$pph4_thresh==pph4_threshold,], pph4_threshold)
output_file <- paste0(visualize_coloc_dir, "number_of_colocalizations_stratefied_by_component_and_type_pph4_thresh_", pph4_threshold, "_bar_plot.pdf")
ggsave(joint_bar_plot, file=output_file, width=7.2, height=9.0, units="in")


