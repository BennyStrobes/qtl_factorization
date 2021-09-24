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


make_single_manhatten_plot <- function(df, study_name, coloring) {
  df$position_real <- as.numeric(df$position)
  plotter <- ggplot(df) + 
             geom_point(aes(x=position_real, y=neg_log_pvalue), size=.1, color=coloring) +
             scale_color_manual(values=c(coloring)) +
             gtex_v8_figure_theme() + 
             labs(x="Position",y="-log10(pvalue)", title=study_name)
  return(plotter)
}

make_single_manhatten_plot_pretty <- function(df, study_name) {
  df$position_real <- as.numeric(df$position)
  plotter <- ggplot(df) + 
             geom_point(aes(x=position_real, y=neg_log_pvalue), size=.1) +
             gtex_v8_figure_theme() + 
             labs(x="",y="-log10(pvalue)", title="") +
             theme(axis.title.x=element_blank()) +
             theme(title = element_blank())
  return(plotter)
}

make_single_upside_down_manhatten_plot <- function(df, study_name) {
  df$position_real <- as.numeric(df$position)
  plotter <- ggplot(df) + 
             geom_point(aes(x=position_real, y=neg_log_pvalue), size=.1) +
             gtex_v8_figure_theme() + 
             labs(x="",y="-log10(pvalue)", title="") +
             scale_y_reverse() + scale_x_continuous(position="top") +
             theme(axis.title.x=element_blank()) +
             theme(axis.text.x=element_blank()) +
             theme(title = element_blank())
  return(plotter)
}

plot_coloc_pretty <- function(eqtl_df, gwas_df, eqtl_study_name, gwas_study_name) {
	eqtl_manhatten <- make_single_upside_down_manhatten_plot(eqtl_df, eqtl_study_name)
	gwas_manhatten <- make_single_manhatten_plot_pretty(gwas_df, gwas_study_name)

	joint <- ggdraw() +
		draw_plot(gwas_manhatten+ theme(legend.position="none"),0,.48,1,.53) +
		draw_plot(eqtl_manhatten + theme(legend.position="none"),0,0,1,.51) 


	return(joint)
}


plot_coloc <- function(eqtl_df, gwas_df, eqtl_study_name, gwas_study_name) {
	eqtl_manhatten <- make_single_manhatten_plot(eqtl_df, eqtl_study_name, "dodgerblue3")
	gwas_manhatten <- make_single_manhatten_plot(gwas_df, gwas_study_name, "chartreuse4")

	joint <- plot_grid(gwas_manhatten, eqtl_manhatten, ncol=1)

	return(joint)
}




output_root <- args[1]
eqtl_study_name <- args[2]
gwas_study_name <- args[3]

test_info_file = paste0(output_root, "test_info.txt")
test_info <- read.table(test_info_file, header=TRUE)

pph4_vec <- c()

num_genes <- dim(test_info)[1]
for (gene_index in 1:num_genes) {
	gene_name <- test_info$gene_name[gene_index]
	chrom_num <- test_info$chrom_num[gene_index]
	eqtl_data_file <- as.character(test_info$eqtl_data_file[gene_index])
	gwas_data_file <- as.character(test_info$gwas_data_file[gene_index])
	eqtl_df <- read.table(eqtl_data_file, header=TRUE, sep="\t")
	gwas_df <- read.table(gwas_data_file, header=TRUE, sep="\t")
	eqtl_df$pvalue <- pnorm( -abs( eqtl_df$beta/sqrt(eqtl_df$varbeta) ) ) * 2
	gwas_df$pvalue <- pnorm( -abs( gwas_df$beta/sqrt(gwas_df$varbeta) ) ) * 2
	eqtl_df$neg_log_pvalue = -log10(eqtl_df$pvalue + 1e-50)
	gwas_df$neg_log_pvalue = -log10(gwas_df$pvalue + 1e-50)

	eqtl_data <- as.list(eqtl_df)
	gwas_data <- as.list(gwas_df)
	eqtl_data$snp = as.character(eqtl_data$snp)
	gwas_data$snp = as.character(gwas_data$snp)

	eqtl_data$type = "quant"
	gwas_data$type = "quant"
	eqtl_data$sdY = 1.0
	gwas_data$sdY = 1.0

	my.res <- coloc.abf(dataset1=eqtl_data, dataset2=gwas_data)

	pph4 <- my.res$summary[6]

	if ((pph4 > .95)) {
		pp <- plot_coloc(eqtl_df, gwas_df, eqtl_study_name, gwas_study_name)
		ggsave(pp, file=paste0(output_root, gene_name, "_viz.pdf"), width=7.2, height=5.0, units="in")

	}

	pph4_vec <- c(pph4_vec, pph4)
}

test_info$pph4 <- pph4_vec
write.table(test_info, file = paste0(output_root, "test_results.txt"), quote = FALSE, sep = "\t",row.names = FALSE)
print(sum(pph4_vec > .7))
print(sum(pph4_vec > .9))
print(sum(pph4_vec > .95))
print(sum(pph4_vec > .99))