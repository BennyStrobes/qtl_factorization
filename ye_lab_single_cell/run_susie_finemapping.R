args = commandArgs(trailingOnly=TRUE)
library(coloc)
library(grid)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
options(bitmapType = 'cairo', device = 'pdf')

gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


make_single_manhatten_plot <- function(df) {
  df$position_real <- as.numeric(df$position)
  plotter <- ggplot(df) + 
             geom_point(aes(x=position_real, y=neg_log_pvalue), size=.1) +
             gtex_v8_figure_theme() + 
             labs(x="Position",y="-log10(pvalue)")
  return(plotter)
}











output_processed_data_root <- args[1]
output_results_root <- args[2]
output_visualization_root <- args[3]


test_info_file = paste0(output_processed_data_root, "test_info.txt")
test_info <- read.table(test_info_file, header=TRUE)


gene_index <- 14

	gene_name <- test_info$gene_name[gene_index]
	chrom_num <- test_info$chrom_num[gene_index]
	eqtl_data_file <- as.character(test_info$eqtl_data_file[gene_index])
	eqtl_df <- read.table(eqtl_data_file, header=TRUE, sep="\t")
	eqtl_df$pvalue <- pnorm( -abs( eqtl_df$beta/sqrt(eqtl_df$varbeta) ) ) * 2
	eqtl_df$neg_log_pvalue = -log10(eqtl_df$pvalue + 1e-50)
	eqtl_ld_file <- as.character(test_info$ld_mat_file[gene_index])
	eqtl_ld <- read.table(eqtl_ld_file, header=FALSE, sep="\t")
	colnames(eqtl_ld) = as.character(eqtl_df$snp)
	rownames(eqtl_ld) = as.character(eqtl_df$snp)

	eqtl_data <- as.list(eqtl_df)
	eqtl_data$snp = as.character(eqtl_data$snp)
	eqtl_data$type = "quant"
	eqtl_data$sdY = 1.0
	eqtl_data$LD = as.matrix(eqtl_ld)
	check_dataset(eqtl_data,req="LD")
	susie_results = runsusie(eqtl_data)
	print(summary(susie_results))


	#pp <- make_single_manhatten_plot(eqtl_df)
	#ggsave(pp, file=paste0(output_visualization_root, gene_name, "_", gene_index, "_viz.pdf"), width=7.2, height=5.0, units="in")
