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


get_ordered_tfs_by_odds_ratio <- function(tf_enrichment_df) {
	unique_tfs <- as.character(unique(tf_enrichment_df$TF))
	ordered_tfs <- as.character(tf_enrichment_df$TF)

	means <- c()
	for (tf_iter in 1:length(unique_tfs)) {
		unique_tf = unique_tfs[tf_iter]
		indices = ordered_tfs == unique_tf 
		means <- c(means, mean(tf_enrichment_df$odds_ratio[indices]))
	}

	indices <- order(means)

	return(unique_tfs[indices])

}

make_tf_odds_ratio_boxplot <- function(df) {
	p <- ggplot(data=df, aes(x=TF, y=odds_ratio)) +
		geom_boxplot(outlier.size=1e-8) +
		gtex_v8_figure_theme() +
		labs(y="odds_ratio", x="TF") +
		geom_hline(yintercept=1) +
		theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
	return(p)
}

make_tf_pvalue_barplot <- function(df) {
	df$neglog10_pvalue = -log10(df$pvalue + 1e-5)
	p <- ggplot(data=df, aes(x=TF, y=neglog10_pvalue)) +
		geom_bar(stat="identity", color="black", position=position_dodge(), width=.75) +
		gtex_v8_figure_theme() +
		labs(y="-log10(pvalue)", x="Transcription factor") +
		geom_hline(yintercept=-log10(.05/1000), color='red') +
		theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

	return(p)
}


working_dir <- args[1]



enrichments_file <- paste0(working_dir, "remap_enrichment_probabilities.txt")

tf_enrichment_df <- read.table(enrichments_file, header=TRUE)
indices <- order(tf_enrichment_df$pvalue)
tf_enrichment_df$TF = factor(tf_enrichment_df$TF, levels=as.character(tf_enrichment_df$TF[indices]))

output_file <- paste0(working_dir, "pvalue_barplot.pdf")
pvalue_barplot<- make_tf_pvalue_barplot(tf_enrichment_df)
ggsave(pvalue_barplot, file=output_file, width=11.2, height=6.0, units="in")


#output_file <- paste0(working_dir, "odds_ratio_boxplot.pdf")
#orat_boxplot <- make_tf_odds_ratio_boxplot(tf_enrichment_df)
#ggsave(orat_boxplot, file=output_file, width=11.2, height=6.0, units="in")


