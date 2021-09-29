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

neg_log10_pvalue_scatter_comparison <- function(marginal_pvalue, joint_pvalue, factor_num) {

	df <- data.frame(marginal=-log10(marginal_pvalue + 1e-50), joint=-log10(joint_pvalue+1e-50))

	print(cor(df$marginal, df$joint))

	scatter <- ggplot(df, aes(x=marginal, y=joint)) +
  				geom_point(size=.1) +
  				gtex_v8_figure_theme() + 
  				geom_abline() +
  				labs(x="Marginal interaction -log10(pvalue)",y="joint interaction -log10(pvalue)", title=paste0("latent context ", factor_num)) 
  	return(scatter)
}




working_dir <- args[1]


input_file <- paste0(working_dir, "interaction_eqtl_results_0_200_results.txt")
input_data <- read.table(input_file, header=FALSE)


for (factor_num in 1:10) {
marginal_pvalue_column <- 5 + 3*(factor_num-1)
joint_pvalue_column <- 5 + 30 + 3*(factor_num-1)
output_file <- paste0(working_dir, "scatter_of_neg_log_10_pvalues_", factor_num, ".pdf")
scatter <- neg_log10_pvalue_scatter_comparison(input_data[,marginal_pvalue_column], input_data[, joint_pvalue_column], factor_num)
ggsave(scatter, file=output_file, width=7.2, height=5, units="in")

}


