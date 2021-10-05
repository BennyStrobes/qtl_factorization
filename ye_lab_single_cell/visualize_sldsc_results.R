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


extract_study_names <- function(processed_gwas_studies_file) {
	df <- read.table(processed_gwas_studies_file, header=FALSE)
	return(as.character(df$V1))
}

extract_surge_tau_df <- function(study_names, sldsc_results_dir, suffix) {
	study_arr <- c()
	category_arr <- c()
	tau_arr <- c()
	for (study_iter in 1:length(study_names)) {
		study_name <- study_names[study_iter]
		cat_levels <- c()
		# Add surge eqtls
		for (latent_factor_num in 1:10) {
			study_file <- paste0(sldsc_results_dir, study_name, "_surge_", latent_factor_num, suffix)
			study_df <- read.table(study_file, header=TRUE)
			num_rows <- dim(study_df)[1]
			study_df_small <- study_df[98:num_rows,]

			study_arr <- c(study_arr, study_name)
			category_arr <- c(category_arr, paste0("surge_eqtl_", latent_factor_num))
			tau_arr <- c(tau_arr, study_df_small$Coefficient_z.score[1])

			cat_levels <- c(cat_levels, paste0("surge_eqtl_", latent_factor_num))
		}
		
		# Add standard eqtl
		study_file <- paste0(sldsc_results_dir, study_name, "_standard", suffix)
		study_df <- read.table(study_file, header=TRUE)
		num_rows <- dim(study_df)[1]
		study_df_small <- study_df[98:num_rows,]
		study_arr <- c(study_arr, study_name)
		category_arr <- c(category_arr, "standard_eqtl")
		tau_arr <- c(tau_arr, study_df_small$Coefficient_z.score[1])

		cat_levels <- c(cat_levels, "standard_eqtl")
	}

	study_levels_non_blood_immune <- c("Alzheimer", "Bipolar", "CAD", "Schizophrenia", "ukbb_bmi", "ukbb_height", "ukbb_T2D")
	study_levels_immune <- c("Celiac", "Crohns", "IBD", "Lupus", "Multiple_sclerosis", "PBC", "Rheumatoid_Arthritis", "ukbb_eczema", "Ulcerative_Colitis")
	study_levels_blood <- c("ukbb_blood_eosinophil_count", "ukbb_blood_high_light_scatter_reticulotye_count", "ukbb_blood_lymphocyte_count", "ukbb_blood_mean_corpuscular_hemoglobin", "ukbb_blood_monocyte_count", "ukbb_blood_platelet_count", "ukbb_blood_platelet_vol", "ukbb_blood_red_count", "ukbb_blood_white_count")
	study_levels <- c(study_levels_non_blood_immune, study_levels_immune, study_levels_blood)
	
	df <- data.frame(study=factor(study_arr, levels=study_levels), category=factor(category_arr, levels=cat_levels), tau=tau_arr)
	df$study = revalue(df$study, c("ukbb_bmi"="bmi", "ukbb_eczema"="eczema", "ukbb_blood_eosinophil_count"="eosinophil count", "ukbb_blood_high_light_scatter_reticulotye_count"="reticulocyte count", "ukbb_blood_lymphocyte_count"="lymphocyte count", "ukbb_blood_mean_corpuscular_hemoglobin"="corpuscular hemoglobin", "ukbb_blood_monocyte_count"="monocyte count", "ukbb_blood_platelet_count"="platelet count", "ukbb_blood_red_count"="blood red count", "ukbb_blood_white_count"="blood white count", "ukbb_blood_platelet_vol"="blood platelet vol"))

	return(df)

}

extract_surge_enrichment_df <- function(study_names, sldsc_results_dir, suffix) {

	study_arr <- c()
	category_arr <- c()
	prop_snps_arr <- c()
	prop_h2_arr <- c()
	prop_h2_std_err_arr <- c()
	enrichment_arr <- c()
	enrichment_std_err_arr <- c()
	enrichment_p_arr <- c()
	tau_arr = c()


	for (study_iter in 1:length(study_names)) {

		study_name <- study_names[study_iter]
		study_file <- paste0(sldsc_results_dir, study_name, suffix)

		study_df <- read.table(study_file, header=TRUE)
		num_rows <- dim(study_df)[1]
		#study_df_small <- study_df[54:num_rows,]
		study_df_small <- study_df[98:num_rows,]
		new_category <- c()
		old_category <- as.character(study_df_small$Category)
		for (cat_iter in 1:length(old_category)) {
			old_cat <- old_category[cat_iter]
			new_word = strsplit(old_cat, "L2")[[1]][1]
			new_category <- c(new_category, new_word)
		}
		new_category <- as.character(new_category)

		study_arr <- c(study_arr, rep(study_name, length(new_category)))
		category_arr <- c(category_arr, as.character(new_category))
		cat_levels <- as.character(new_category)
		prop_snps_arr <- c(prop_snps_arr, study_df_small$Prop._SNPs)
		prop_h2_arr <- c(prop_h2_arr, study_df_small$Prop._h2)
		prop_h2_std_err_arr <- c(prop_h2_std_err_arr, study_df_small$Prop._h2_std_error)
		enrichment_arr <- c(enrichment_arr, study_df_small$Enrichment)
		enrichment_std_err_arr <- c(enrichment_std_err_arr, study_df_small$Enrichment_std_error)
		enrichment_p_arr <- c(enrichment_p_arr, study_df_small$Enrichment_p)
		tau_arr <- c(tau_arr, study_df_small$Coefficient_z.score)
	}

	study_levels_non_blood_immune <- c("Alzheimer", "Bipolar", "CAD", "Schizophrenia", "ukbb_bmi", "ukbb_height", "ukbb_T2D")
	study_levels_immune <- c("Celiac", "Crohns", "IBD", "Lupus", "Multiple_sclerosis", "PBC", "Rheumatoid_Arthritis", "ukbb_eczema", "Ulcerative_Colitis")
	study_levels_blood <- c("ukbb_blood_eosinophil_count", "ukbb_blood_high_light_scatter_reticulotye_count", "ukbb_blood_lymphocyte_count", "ukbb_blood_mean_corpuscular_hemoglobin", "ukbb_blood_monocyte_count", "ukbb_blood_platelet_count", "ukbb_blood_platelet_vol", "ukbb_blood_red_count", "ukbb_blood_white_count")
	study_levels <- c(study_levels_non_blood_immune, study_levels_immune, study_levels_blood)

	df <- data.frame(study=factor(study_arr, levels=study_levels), category=factor(category_arr, levels=cat_levels), prop_snps=prop_snps_arr, prop_h2=prop_h2_arr, tau=tau_arr, prop_h2_std_err=prop_h2_std_err_arr, enrichment=enrichment_arr, enrichment_std_err=enrichment_std_err_arr, enrichment_p=enrichment_p_arr)

	df$study = revalue(df$study, c("ukbb_bmi"="bmi", "ukbb_eczema"="eczema", "ukbb_blood_eosinophil_count"="eosinophil count", "ukbb_blood_high_light_scatter_reticulotye_count"="reticulocyte count", "ukbb_blood_lymphocyte_count"="lymphocyte count", "ukbb_blood_mean_corpuscular_hemoglobin"="corpuscular hemoglobin", "ukbb_blood_monocyte_count"="monocyte count", "ukbb_blood_platelet_count"="platelet count", "ukbb_blood_red_count"="blood red count", "ukbb_blood_white_count"="blood white count", "ukbb_blood_platelet_vol"="blood platelet vol"))
	return(df)
}


make_enrichment_barplot <- function(df) {
	p <- ggplot(data=df, aes(x=study, y=enrichment, fill=category)) +
		geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
		gtex_v8_figure_theme() +
		theme(legend.position="top") +
		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10)) +
		labs(y="S-LDSC enrichment", fill="") +
		geom_hline(yintercept=1) +
		geom_errorbar(aes(ymin=enrichment-enrichment_std_err, ymax=enrichment+enrichment_std_err), position = position_dodge(), width = .75, size=.2)

	num_categories = length(unique(df$category))
	if (num_categories == 2) {
		p <- p + scale_fill_manual(values=c("hotpink3", "gray53"))
	}
	return(p)
}

make_neglog10_pvalue_barplot <- function(df) {
	df$neglog10_pvalue = -log10(df$enrichment_p)
	p <- ggplot(data=df, aes(x=study, y=neglog10_pvalue, fill=category)) +
		geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
		gtex_v8_figure_theme() +
		theme(legend.position="top") + 
		labs(y="S-LDSC -log10(pvalue)", fill="") +
		geom_hline(yintercept=-log10(.05)) +
		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))
	num_categories = length(unique(df$category))
	if (num_categories == 2) {
		p <- p + scale_fill_manual(values=c("hotpink3", "gray53"))
	}
	return(p)
}

make_cross_trait_neglog10_pvalue_boxplot <- function(df) {
	df$neglog10_pvalue = -log10(df$enrichment_p)
	study_levels_non_blood_immune <- c("Alzheimer", "Bipolar", "CAD", "Schizophrenia", "bmi", "height", "T2D")
	study_levels_immune <- c("Celiac", "Crohns", "IBD", "Lupus", "Multiple_sclerosis", "PBC", "Rheumatoid_Arthritis", "eczema", "Ulcerative_Colitis")
	study_levels_blood <- c("eosinophil_count", "reticulocyte count", "lymphocyte count", "corpuscular hemoglobin", "monocyte count", "platelet count", "blood platelet vol", "blood red count", "blood white count")

	pvalues_arr <- c()
	categories_arr <- c()
	study_groups_arr <- c()

	non_blood_studies = df$study %in% study_levels_non_blood_immune
	pvalues_arr <- c(pvalues_arr, df$neglog10_pvalue[non_blood_studies])
	categories_arr <- c(categories_arr, as.character(df$category[non_blood_studies]))
	study_groups_arr <- c(study_groups_arr, rep("non-blood-immune", sum(non_blood_studies)))

	blood_studies = df$study %in% study_levels_blood
	pvalues_arr <- c(pvalues_arr, df$neglog10_pvalue[blood_studies])
	categories_arr <- c(categories_arr, as.character(df$category[blood_studies]))
	study_groups_arr <- c(study_groups_arr, rep("blood", sum(blood_studies)))

	immune_studies = df$study %in% study_levels_immune
	pvalues_arr <- c(pvalues_arr, df$neglog10_pvalue[immune_studies])
	categories_arr <- c(categories_arr, as.character(df$category[immune_studies]))
	study_groups_arr <- c(study_groups_arr, rep("immune", sum(immune_studies)))

	new_df <- data.frame(pvalues=pvalues_arr, study_group=factor(study_groups_arr), category=factor(categories_arr))
	new_df$category = factor(new_df$category, levels=c("surge_eqtl_1", "surge_eqtl_2", "surge_eqtl_3", "surge_eqtl_4", "surge_eqtl_5", "surge_eqtl_6", "surge_eqtl_7", "surge_eqtl_8", "surge_eqtl_9", "surge_eqtl_10", "standard_eqtl"))
	p <- ggplot(data=new_df, aes(x=category, y=pvalues, fill=study_group)) +
		geom_boxplot() +
		gtex_v8_figure_theme() +
		theme(legend.position="top") + 
		labs(y="S-LDSC -log10(pvalue)", fill="trait category") +
		geom_hline(yintercept=-log10(.05)) +
		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))
	return(p)

}

make_cross_trait_enrichment_boxplot <- function(df) {
	study_levels_non_blood_immune <- c("Alzheimer", "Bipolar", "CAD", "Schizophrenia", "bmi", "height", "T2D")
	study_levels_immune <- c("Celiac", "Crohns", "IBD", "Lupus", "Multiple_sclerosis", "PBC", "Rheumatoid_Arthritis", "eczema", "Ulcerative_Colitis")
	study_levels_blood <- c("eosinophil_count", "reticulocyte count", "lymphocyte count", "corpuscular hemoglobin", "monocyte count", "platelet count", "blood platelet vol", "blood red count", "blood white count")

	pvalues_arr <- c()
	categories_arr <- c()
	study_groups_arr <- c()

	non_blood_studies = df$study %in% study_levels_non_blood_immune
	pvalues_arr <- c(pvalues_arr, df$enrichment[non_blood_studies])
	categories_arr <- c(categories_arr, as.character(df$category[non_blood_studies]))
	study_groups_arr <- c(study_groups_arr, rep("non-blood-immune", sum(non_blood_studies)))

	blood_studies = df$study %in% study_levels_blood
	pvalues_arr <- c(pvalues_arr, df$enrichment[blood_studies])
	categories_arr <- c(categories_arr, as.character(df$category[blood_studies]))
	study_groups_arr <- c(study_groups_arr, rep("blood", sum(blood_studies)))

	immune_studies = df$study %in% study_levels_immune
	pvalues_arr <- c(pvalues_arr, df$enrichment[immune_studies])
	categories_arr <- c(categories_arr, as.character(df$category[immune_studies]))
	study_groups_arr <- c(study_groups_arr, rep("immune", sum(immune_studies)))

	new_df <- data.frame(enrichment=pvalues_arr, study_group=factor(study_groups_arr), category=factor(categories_arr))
	new_df$enrichment[new_df$enrichment < 0.0] = new_df$enrichment[new_df$enrichment < 0.0]*0.0
	new_df$category = factor(new_df$category, levels=c("surge_eqtl_1", "surge_eqtl_2", "surge_eqtl_3", "surge_eqtl_4", "surge_eqtl_5", "surge_eqtl_6", "surge_eqtl_7", "surge_eqtl_8", "surge_eqtl_9", "surge_eqtl_10", "standard_eqtl"))
	p <- ggplot(data=new_df, aes(x=category, y=enrichment, fill=study_group)) +
		geom_boxplot() +
		gtex_v8_figure_theme() +
		theme(legend.position="top") + 
		labs(y="S-LDSC Enrichment", fill="trait category") +
		geom_hline(yintercept=1) +
		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))
	return(p)

}

make_cross_trait_tau_boxplot <- function(df) {
	study_levels_non_blood_immune <- c("Alzheimer", "Bipolar", "CAD", "Schizophrenia", "bmi", "height", "T2D")
	study_levels_immune <- c("Celiac", "Crohns", "IBD", "Lupus", "Multiple_sclerosis", "PBC", "Rheumatoid_Arthritis", "eczema", "Ulcerative_Colitis")
	study_levels_blood <- c("eosinophil_count", "reticulocyte count", "lymphocyte count", "corpuscular hemoglobin", "monocyte count", "platelet count", "blood platelet vol", "blood red count", "blood white count")

	pvalues_arr <- c()
	categories_arr <- c()
	study_groups_arr <- c()

	non_blood_studies = df$study %in% study_levels_non_blood_immune
	pvalues_arr <- c(pvalues_arr, df$tau[non_blood_studies])
	categories_arr <- c(categories_arr, as.character(df$category[non_blood_studies]))
	study_groups_arr <- c(study_groups_arr, rep("non-blood-immune", sum(non_blood_studies)))

	blood_studies = df$study %in% study_levels_blood
	pvalues_arr <- c(pvalues_arr, df$tau[blood_studies])
	categories_arr <- c(categories_arr, as.character(df$category[blood_studies]))
	study_groups_arr <- c(study_groups_arr, rep("blood", sum(blood_studies)))

	immune_studies = df$study %in% study_levels_immune
	pvalues_arr <- c(pvalues_arr, df$tau[immune_studies])
	categories_arr <- c(categories_arr, as.character(df$category[immune_studies]))
	study_groups_arr <- c(study_groups_arr, rep("immune", sum(immune_studies)))

	new_df <- data.frame(tau=pvalues_arr, study_group=factor(study_groups_arr), category=factor(categories_arr))
	new_df$category = factor(new_df$category, levels=c("surge_eqtl_1", "surge_eqtl_2", "surge_eqtl_3", "surge_eqtl_4", "surge_eqtl_5", "surge_eqtl_6", "surge_eqtl_7", "surge_eqtl_8", "surge_eqtl_9", "surge_eqtl_10", "standard_eqtl"))
	p <- ggplot(data=new_df, aes(x=category, y=tau, fill=study_group)) +
		geom_boxplot() +
		gtex_v8_figure_theme() +
		theme(legend.position="top") + 
		labs(y="S-LDSC Tau", fill="trait category") +
		geom_hline(yintercept=0) +
		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))
	return(p)

}

make_cross_trait_enrichment_violinplot <- function(df) {
	study_levels_non_blood_immune <- c("Alzheimer", "Bipolar", "CAD", "Schizophrenia", "bmi", "height", "T2D")
	study_levels_immune <- c("Celiac", "Crohns", "IBD", "Lupus", "Multiple_sclerosis", "PBC", "Rheumatoid_Arthritis", "eczema", "Ulcerative_Colitis")
	study_levels_blood <- c("eosinophil_count", "reticulocyte count", "lymphocyte count", "corpuscular hemoglobin", "monocyte count", "platelet count", "blood platelet vol", "blood red count", "blood white count")

	pvalues_arr <- c()
	categories_arr <- c()
	study_groups_arr <- c()

	non_blood_studies = df$study %in% study_levels_non_blood_immune
	pvalues_arr <- c(pvalues_arr, df$enrichment[non_blood_studies])
	categories_arr <- c(categories_arr, as.character(df$category[non_blood_studies]))
	study_groups_arr <- c(study_groups_arr, rep("non-blood-immune", sum(non_blood_studies)))

	blood_studies = df$study %in% study_levels_blood
	pvalues_arr <- c(pvalues_arr, df$enrichment[blood_studies])
	categories_arr <- c(categories_arr, as.character(df$category[blood_studies]))
	study_groups_arr <- c(study_groups_arr, rep("blood", sum(blood_studies)))

	immune_studies = df$study %in% study_levels_immune
	pvalues_arr <- c(pvalues_arr, df$enrichment[immune_studies])
	categories_arr <- c(categories_arr, as.character(df$category[immune_studies]))
	study_groups_arr <- c(study_groups_arr, rep("immune", sum(immune_studies)))

	new_df <- data.frame(enrichment=pvalues_arr, study_group=factor(study_groups_arr), category=factor(categories_arr))
	new_df$enrichment[new_df$enrichment < 0.0] = new_df$enrichment[new_df$enrichment < 0.0]*0.0
	new_df$category = factor(new_df$category, levels=c("surge_eqtl_1", "surge_eqtl_2", "surge_eqtl_3", "surge_eqtl_4", "surge_eqtl_5", "surge_eqtl_6", "surge_eqtl_7", "surge_eqtl_8", "surge_eqtl_9", "surge_eqtl_10", "standard_eqtl"))
	p <- ggplot(data=new_df, aes(x=category, y=enrichment, fill=study_group)) +
		geom_violin() +
		gtex_v8_figure_theme() +
		theme(legend.position="top") + 
		labs(y="S-LDSC Enrichment", fill="trait category") +
		geom_hline(yintercept=1) +
		theme(axis.text.x = element_text(angle = 90,hjust=1, vjust=.5, size=10))
	return(p)

}







processed_gwas_studies_file = args[1]
sldsc_results_dir = args[2]
sldsc_visualization_dir = args[3]


# Extract gwas study names
study_names <- extract_study_names(processed_gwas_studies_file)

# Extract data frame of enrichment
sldsc_surge_enrichment_df = extract_surge_enrichment_df(study_names, sldsc_results_dir, "_surge_eqtls_1e-05_and_baselineLD.results")
sldsc_surge_tau_df = extract_surge_tau_df(study_names, sldsc_results_dir, "_eqtls_1e-05_only_and_baselineLD.results")
#sldsc_surge_enrichment_df = extract_surge_enrichment_df(study_names, sldsc_results_dir, "_surge_eqtls_1e-05_and_baseline.results")
#sldsc_joint_surge_enrichment_df = extract_surge_enrichment_df(study_names, sldsc_results_dir, "_joint_surge_egenes_05_and_baseline.results")


##############################
# Make enrichment bar plot
##############################
output_file <- paste0(sldsc_visualization_dir, "surge_eqtls_1e-05_and_baselineLD_enrichment_barplot.pdf")
enrichment_barplot <- make_enrichment_barplot(sldsc_surge_enrichment_df)
ggsave(enrichment_barplot, file=output_file, width=12.2, height=6.0, units="in")



##############################
# Make -log10 pvalue bar plot
##############################
output_file <- paste0(sldsc_visualization_dir, "surge_eqtls_1e-05_and_baselineLD_neglog10_pvalue_barplot.pdf")
neglog10_pvalue_barplot <- make_neglog10_pvalue_barplot(sldsc_surge_enrichment_df)
ggsave(neglog10_pvalue_barplot, file=output_file, width=12.2, height=6.0, units="in")




##############################
# Make cross-trait -log10 pvalue box plot
##############################
output_file <- paste0(sldsc_visualization_dir, "surge_eqtls_1e-05_and_baselineLD_cross_trait_neglog10_pvalue_boxplot.pdf")
x_trait_neglog10_pvalue_boxplot <- make_cross_trait_neglog10_pvalue_boxplot(sldsc_surge_enrichment_df)
ggsave(x_trait_neglog10_pvalue_boxplot, file=output_file, width=12.2, height=6.0, units="in")

##############################
# Make cross-trait enrichment box plot
##############################
output_file <- paste0(sldsc_visualization_dir, "surge_eqtls_1e-05_and_baselineLD_cross_trait_enrichment_boxplot.pdf")
x_trait_enrichment_boxplot <- make_cross_trait_enrichment_boxplot(sldsc_surge_enrichment_df)
ggsave(x_trait_enrichment_boxplot, file=output_file, width=12.2, height=6.0, units="in")


##############################
# Make cross-trait enrichment violin plot
##############################
output_file <- paste0(sldsc_visualization_dir, "surge_eqtls_1e-05_and_baselineLD_cross_trait_enrichment_violinplot.pdf")
x_trait_enrichment_boxplot <- make_cross_trait_enrichment_violinplot(sldsc_surge_enrichment_df)
ggsave(x_trait_enrichment_boxplot, file=output_file, width=12.2, height=6.0, units="in")



##############################
# Make cross-trait tau box plot
##############################
output_file <- paste0(sldsc_visualization_dir, "surge_eqtls_1e-05_one_at_a_time_and_baselineLD_cross_trait_tau_boxplot.pdf")
x_trait_tau_boxplot <- make_cross_trait_tau_boxplot(sldsc_surge_tau_df)
ggsave(x_trait_tau_boxplot, file=output_file, width=12.2, height=6.0, units="in")

output_file <- paste0(sldsc_visualization_dir, "surge_eqtls_1e-05_and_baselineLD_cross_trait_tau_boxplot.pdf")
x_trait_tau_boxplot <- make_cross_trait_tau_boxplot(sldsc_surge_enrichment_df)
ggsave(x_trait_tau_boxplot, file=output_file, width=12.2, height=6.0, units="in")

