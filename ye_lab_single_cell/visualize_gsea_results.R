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

make_single_component_odds_ratio_plot <- function(enrichment_df_subset, component_num) {
	enrichment_df_subset$adjusted_log_pvalue = -log10(enrichment_df_subset$adjusted_pvalue)
	#enrichment_df$log_odds_ratio = log(enrichment_df$odds_ratio)
	enrichment_df_subset$pathway=factor(enrichment_df_subset$pathway, levels=sort(enrichment_df_subset$pathway, decreasing = TRUE))
	p <- ggplot() + 
					geom_point(data=enrichment_df_subset, mapping=aes(y=pathway, x=odds_ratio)) +
					labs(y = "Pathway", x = "odds ratio", title=paste0("Latent context ", component_num)) +
					gtex_v8_figure_theme() + theme(axis.text.y = element_text(size=8))
	return(p)
}

make_cross_component_odds_ratio_plot <- function(enrichment_df, width) {
	component_num = 3
	p3 <- make_single_component_odds_ratio_plot(enrichment_df[enrichment_df$component_num == component_num,], component_num+1)
	component_num = 4
	p4 <- make_single_component_odds_ratio_plot(enrichment_df[enrichment_df$component_num == component_num,], component_num+1)
	component_num = 5
	p5 <- make_single_component_odds_ratio_plot(enrichment_df[enrichment_df$component_num == component_num,], component_num+1)
	component_num = 6
	p6 <- make_single_component_odds_ratio_plot(enrichment_df[enrichment_df$component_num == component_num,], component_num+1)
	component_num = 7
	p7 <- make_single_component_odds_ratio_plot(enrichment_df[enrichment_df$component_num == component_num,], component_num+1)
	component_num = 8
	p8 <- make_single_component_odds_ratio_plot(enrichment_df[enrichment_df$component_num == component_num,], component_num+1)
	component_num = 9
	p9 <- make_single_component_odds_ratio_plot(enrichment_df[enrichment_df$component_num == component_num,], component_num+1)

  	merged = plot_grid(p3 + theme(legend.position='none'), p4 + theme(legend.position='none', axis.text.y=element_blank()) +theme(axis.title.y=element_blank()), p5 + theme(legend.position='none', axis.text.y=element_blank()) +theme(axis.title.y=element_blank()), p6 + theme(legend.position='none', axis.text.y=element_blank()) +theme(axis.title.y=element_blank()), p7 + theme(legend.position='none', axis.text.y=element_blank()) +theme(axis.title.y=element_blank()), p8 + theme(legend.position='none', axis.text.y=element_blank()) +theme(axis.title.y=element_blank()), p9 + theme(legend.position='none', axis.text.y=element_blank()) +theme(axis.title.y=element_blank()), nrow=1, rel_widths=c(width,1,1,1,1,1,1))
  	return(merged)
}




gene_set_enrichment_output_stem <- args[1]




gsea_type <- "enrichments_c5_bp"

enrichment_file <- paste0(gene_set_enrichment_output_stem, "cross_component_significant_enrichments_", gsea_type, ".txt")
enrichment_df <- read.table(enrichment_file, sep='\t', header=TRUE)


output_file <- paste0(gene_set_enrichment_output_stem, "cross_component_odds_ratio_plot_", gsea_type, ".pdf")
odds_ratio_plot <- make_cross_component_odds_ratio_plot(enrichment_df, 4)
ggsave(odds_ratio_plot, file=output_file, width=14.2, height=6.0, units="in")


gsea_type <- "enrichments_c2_cp_kegg"
enrichment_file <- paste0(gene_set_enrichment_output_stem, "cross_component_significant_enrichments_", gsea_type, ".txt")
enrichment_df <- read.table(enrichment_file, sep='\t', header=TRUE)


output_file <- paste0(gene_set_enrichment_output_stem, "cross_component_odds_ratio_plot_", gsea_type, ".pdf")
odds_ratio_plot <- make_cross_component_odds_ratio_plot(enrichment_df, 3.5)
ggsave(odds_ratio_plot, file=output_file, width=14.2, height=6.0, units="in")
