args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(cowplot)
library(umap)
library(ggplot2)
library(RColorBrewer)
options(bitmapType = 'cairo', device = 'pdf')

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}



make_neg_log10_pvalue_scatter <- function(pvalue1, pvalue2, label1, label2) {
	neg_log_pvalue1 = -log10(pvalue1 + 1e-50)
	neg_log_pvalue2 = -log10(pvalue2 + 1e-50)
	corry = cor(neg_log_pvalue1, neg_log_pvalue2)
	df <- data.frame(neg_log_pvalue1=neg_log_pvalue1, neg_log_pvalue2=neg_log_pvalue2)
	p <- ggplot(df, aes(x=neg_log_pvalue1, y=neg_log_pvalue2)) + geom_point(alpha=.2) +
		labs(x=label1, y=label2, title=paste0("Pearson Correlation: ", corry)) + 
		geom_abline(color="blue") + 
		figure_theme()
	return(p)
}


make_paired_neg_log10_pvalue_qq_plots <- function(pvalue1, pvalue2, label1, label2) {
	neg_log_pvalue1 = -log10(pvalue1 + 1e-50)
	neg_log_pvalue2 = -log10(pvalue2 + 1e-50)
	uniform_pvalues = runif(length(neg_log_pvalue1))
	neg_log_uniform_pvalues = -log10(uniform_pvalues + 1e-50)

	real_pvalues <- c(sort(neg_log_pvalue1), sort(neg_log_pvalue2))
	null_pvalues <- c(sort(neg_log_uniform_pvalues), sort(neg_log_uniform_pvalues))
	pvalue_types <- c(rep(label1, length(neg_log_pvalue1)), rep(label2, length(neg_log_pvalue2)))

	df <- data.frame(real_pvalues=real_pvalues, null_pvalues=null_pvalues, pvalue_types=factor(pvalue_types))

	p <- ggplot(df, aes(x=null_pvalues, y=real_pvalues, color=pvalue_types)) + geom_point() +
		labs(x="Null", y="Real", color="Interaction variable") + 
		geom_abline(color="blue") + 
		figure_theme()
	return(p)
}

make_num_hits_per_latent_factor_barplot <- function(num_hits_per_latent_factor, pvalue_threshold) {
	num_hits_per_latent_factor <- num_hits_per_latent_factor[num_hits_per_latent_factor$nominal_pvalue_threshold==pvalue_threshold,]
	num_hits_per_latent_factor$latent_factor = factor(num_hits_per_latent_factor$latent_factor)
	p<-ggplot(data=num_hits_per_latent_factor, aes(x=latent_factor, y=number_of_hits)) +
  		geom_bar(stat="identity") + 
  		figure_theme() +
  		labs(title=paste0("Nominal coefficient pvalue threshold: ", pvalue_threshold), x="Latent factor", y="Number of hits") +
  		theme(axis.text.x = element_text(size=8,angle = 90, vjust = 0.5, hjust=1))
  	return(p)
}

make_num_hits_per_latent_factor_barplot_shell <- function(num_hits_per_latent_factor) {
	p1 <- make_num_hits_per_latent_factor_barplot(num_hits_per_latent_factor, 1e-4)
	p2 <- make_num_hits_per_latent_factor_barplot(num_hits_per_latent_factor, 1e-6)
	p3 <- make_num_hits_per_latent_factor_barplot(num_hits_per_latent_factor, 1e-8)
	combined <- plot_grid(p1, p2, p3, ncol=1)
	return(combined)
}

make_num_hits_per_tissue_barplot <- function(num_hits_per_latent_factor, pvalue_threshold, tissue_names) {
	num_hits_per_latent_factor <- num_hits_per_latent_factor[num_hits_per_latent_factor$nominal_pvalue_threshold==pvalue_threshold,]
	num_hits_per_latent_factor$latent_factor = factor(num_hits_per_latent_factor$latent_factor)
	p<-ggplot(data=num_hits_per_latent_factor, aes(x=latent_factor, y=number_of_hits)) +
  		geom_bar(stat="identity") + 
  		figure_theme() +
  		labs(title=paste0("Nominal coefficient pvalue threshold: ", pvalue_threshold), x="Tissue", y="Number of hits") +
  		scale_x_discrete(labels=as.character(tissue_names)) +
  		theme(axis.text.x = element_text(size=8,angle = 90, vjust = 0.5, hjust=1))
  	return(p)
}


make_num_hits_per_tissue_barplot_shell <- function(num_hits_per_latent_factor, tissue_names) {
	p1 <- make_num_hits_per_tissue_barplot(num_hits_per_latent_factor, 1e-4, tissue_names)
	p2 <- make_num_hits_per_tissue_barplot(num_hits_per_latent_factor, 1e-6, tissue_names)
	p3 <- make_num_hits_per_tissue_barplot(num_hits_per_latent_factor, 1e-8, tissue_names)
	combined <- plot_grid(p1, p2, p3, ncol=1)
	return(combined)
}

make_num_hits_per_latent_factor_quantile_stratefied_barplot <- function(num_hits_per_latent_factor, quantile_bin) {
	num_hits_per_latent_factor <- num_hits_per_latent_factor[num_hits_per_latent_factor$quantile_bin==quantile_bin,]
	num_hits_per_latent_factor$latent_factor = factor(num_hits_per_latent_factor$latent_factor)
	p<-ggplot(data=num_hits_per_latent_factor, aes(x=latent_factor, y=number_of_hits)) +
  		geom_bar(stat="identity") + 
  		figure_theme() +
  		labs(title=paste0("Quantile bin: ", quantile_bin), x="Latent factor", y="Number of hits") +
  		theme(axis.text.x = element_text(size=8,angle = 90, vjust = 0.5, hjust=1))
  	return(p)
}

make_num_hits_per_latent_factor_quantile_stratefied_barplot_shell <- function(num_hits_per_latent_factor) {
	p1 <- make_num_hits_per_latent_factor_quantile_stratefied_barplot(num_hits_per_latent_factor, 0)
	p2 <- make_num_hits_per_latent_factor_quantile_stratefied_barplot(num_hits_per_latent_factor, 1)
	p3 <- make_num_hits_per_latent_factor_quantile_stratefied_barplot(num_hits_per_latent_factor, 2)
	p4 <- make_num_hits_per_latent_factor_quantile_stratefied_barplot(num_hits_per_latent_factor, 3)
	p5 <- make_num_hits_per_latent_factor_quantile_stratefied_barplot(num_hits_per_latent_factor, 4)
	combined <- plot_grid(p1, p2, p3, p4, p5, ncol=1)
	return(combined)
}


working_dir <- args[1]
tissues_file <- args[2]
visualization_dir <- args[3]



###################
# Input data
###################

# Subsset of results from latent-factor interaction analysis
lf_interaction_results_file <- paste0(working_dir, "cross_tissue_latent_factor_interaction_eqtl_results_merged.txt")
lf_interaction_pvalues <- read.table(lf_interaction_results_file, header=FALSE, sep="\t", nrows=50000)$V3

# Subsset of results from known-tissue interaction analysis
kt_interaction_results_file <- paste0(working_dir, "cross_tissue_known_tissue_interaction_eqtl_results_merged.txt")
kt_interaction_pvalues <- read.table(kt_interaction_results_file, header=FALSE, sep="\t", nrows=50000)$V3

# File containing tally of number of FDR hits driven by each latent factor
num_hits_per_latent_factor_file <- paste0(working_dir, "cross_tissue_latent_factor_interaction_eqtl_results_number_of_hits_per_latent_factor.txt")
num_hits_per_latent_factor <- read.table(num_hits_per_latent_factor_file, header=TRUE, sep="\t")

# File containing tally of number of FDR hits driven by each tissue
num_hits_per_tissue_file <- paste0(working_dir, "cross_tissue_known_tissue_interaction_eqtl_results_number_of_hits_per_latent_factor.txt")
num_hits_per_tissue <- read.table(num_hits_per_tissue_file, header=TRUE, sep="\t")

# File containing tally of number of FDR hits driven by each latent factor stratefied by pvalue quantiles
num_hits_per_latent_factor_quantile_stratefied_file <- paste0(working_dir, "cross_tissue_latent_factor_interaction_eqtl_results_number_of_hits_per_latent_factor_stratefied_by_5_quantiles.txt")
num_hits_per_latent_factor_quantile_stratefied <- read.table(num_hits_per_latent_factor_quantile_stratefied_file, header=TRUE, sep="\t")


# File containing tissue names
tissue_names <- read.table(tissues_file, header=FALSE)$V1


###################
# Make plots
###################

# Make scatter plot of lf interactions vs kt interactions
output_file <- paste0(visualization_dir, "scatterplot_of_lf_and_kt_interaction_pvalues.pdf")
# scatter <- make_neg_log10_pvalue_scatter(lf_interaction_pvalues, kt_interaction_pvalues, "-log10(latent factor interaction pvalue)", "-log10(known tissue interaction pvalue)")
# ggsave(scatter, file=output_file, width=7.2, height=6.0, units="in")


# Make qq plot of lf interactions vs kt interactions
output_file <- paste0(visualization_dir, "qqplot_of_lf_and_kt_interaction_pvalues.pdf")
#qq_plot <- make_paired_neg_log10_pvalue_qq_plots(lf_interaction_pvalues, kt_interaction_pvalues, "latent factor interaction", "known tissue interaction")
#ggsave(qq_plot, file=output_file, width=7.2, height=6.0, units="in")



# Make barplot of number of hits per latent factor
output_file <- paste0(visualization_dir, "barplot_of_num_hits_per_lf.pdf")
#barplot <- make_num_hits_per_latent_factor_barplot_shell(num_hits_per_latent_factor)
#ggsave(barplot, file=output_file, width=9.2, height=5.0, units="in")


# Make barplot of number of hits per latent factor
output_file <- paste0(visualization_dir, "barplot_of_num_hits_per_tissue.pdf")
#barplot <- make_num_hits_per_tissue_barplot_shell(num_hits_per_tissue, tissue_names)
#ggsave(barplot, file=output_file, width=7.2, height=14.0, units="in")


# Make barplot of number of hits per latent factor
output_file <- paste0(visualization_dir, "barplot_of_num_hits_per_lf_quantile_stratefied.pdf")
barplot <- make_num_hits_per_latent_factor_quantile_stratefied_barplot_shell(num_hits_per_latent_factor_quantile_stratefied)
ggsave(barplot, file=output_file, width=9.2, height=8.0, units="in")

