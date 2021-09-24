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


make_enrichment_boxplot_with_columns_for_trait_colored_by_context <- function(df) {

	boxplot <- ggplot(df, aes(x=study, y=odds_ratio, fill=context)) + geom_boxplot(outlier.size = .001) +
				gtex_v8_figure_theme() + 
	        	labs(x="GWAS Study", y = "pi1_real - pi1_background", fill="Surge latent context") +
	        	theme(legend.position="bottom") +
	           	geom_hline(yintercept = 0) + 
	           	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

	return(boxplot)
}









gwas_overlap_dir <- args[1]

#################################################
############# Load in data into nice df
#################################################
latent_contexts <- c(1,2,3,4,6,7)
odds_ratio_vec <- c()
study_name_vec <- c()
latent_context_vec <- c()
pvalue_vec <- c()
for (iter in 1:length(latent_contexts)) {
	latent_context <- latent_contexts[iter]
	file_name <- paste0(gwas_overlap_dir, "gwas_overlap_surge_context_", latent_context, "_cross_study_gwas_enrichment.txt")
	context_data <- read.table(file_name, header=TRUE)
	
	study_names <- colnames(context_data)
	num_studies = length(study_names)
	for (study_iter in 1:num_studies) {
		study_name <- study_names[study_iter]
		orat <- context_data[,study_iter]
		orat[which(!is.finite(orat))] <- 10
		odds_ratio_vec <- c(odds_ratio_vec, orat)
		study_name_vec <- c(study_name_vec, rep(study_name, length(orat)))
		latent_context_vec <- c(latent_context_vec, rep(latent_context, length(orat)))
	}
}
df <- data.frame(study=factor(study_name_vec, levels=c("AD", "BP", "BMI", "CAD", "Height", "SCZ", "T2D", "Celiac", "Crohns", "IBD", "Lupus", "MS", "PBC", "RA", "UC", "immune_average", "non_immune_average")), context=factor(latent_context_vec, levels=latent_contexts), log_odds_ratio=log(odds_ratio_vec), odds_ratio=odds_ratio_vec)



#################################################
############# Start visualizing
#################################################
######################################
# Make enrichment boxplot colored by context number stratefied by disease type
#######################################
output_file <- paste0(gwas_overlap_dir, "enrichment_pvalue_boxplot_with_columns_for_trait_colored_by_context.pdf")
boxplot <- make_enrichment_boxplot_with_columns_for_trait_colored_by_context(df)
ggsave(boxplot, file=output_file, width=13.2, height=4.5, units="in")

output_file <- paste0(gwas_overlap_dir, "enrichment_pvalue_boxplot_with_columns_for_average_trait_colored_by_context.pdf")
boxplot <- make_enrichment_boxplot_with_columns_for_trait_colored_by_context(df[(as.character(df$study) == "immune_average" | as.character(df$study) == "non_immune_average"), ])
ggsave(boxplot, file=output_file, width=7.2, height=4.5, units="in")







