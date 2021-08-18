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


load_in_data <- function(input_data_dir) {
	tau_vec <- c()
	version_vec <- c()
	af_vec <- c()
	seed_vec <- c()
	test_loaded_vec <- c()


	seeds <- c(1, 2, 3, 4, 5)
	afs <- c(".1",".2",".3",".4",".5",".6",".7",".8",".9")
	for (seed_iter in 1:length(seeds)) {
		for (af_iter in 1:length(afs)) {
			seed <- seeds[seed_iter]
			af <- afs[af_iter]
			file_name <- paste0(input_data_dir, "results_af_", af, "_version_heteroskedasticity_seed_", seed, "_analytics.txt")
			aa <- read.table(file_name, header=FALSE, sep='\t')
			tau <- aa[1,1]
			test_loaded <- as.character(aa[1,2])

			tau_vec <- c(tau_vec, tau)
			version_vec <- c(version_vec, "heteroskedasticity")
			af_vec <- c(af_vec, as.numeric(af))
			seed_vec <- c(seed_vec, seed)
			test_loaded_vec <- c(test_loaded_vec, test_loaded)
		}
	}

	seeds <- c(1, 2, 3, 4, 5)
	af <- ".1"
	for (seed_iter in 1:length(seeds)) {
		seed <- seeds[seed_iter]
		file_name <- paste0(input_data_dir, "results_af_", af, "_version_none_seed_", seed, "_analytics.txt")
		aa <- read.table(file_name, header=FALSE, sep='\t')
		tau <- aa[1,1]
		test_loaded <- as.character(aa[1,2])

		tau_vec <- c(tau_vec, tau)
		version_vec <- c(version_vec, "none")
		af_vec <- c(af_vec, as.numeric(af))
		seed_vec <- c(seed_vec, seed)
		test_loaded_vec <- c(test_loaded_vec, test_loaded)
	}


	df <- data.frame(tau=tau_vec, version=version_vec, af=af_vec, seed=seed_vec, test_loaded=factor(test_loaded_vec))
	return(df)
}

make_af_vs_tau_scatter_plot_colored_by_whether_test_was_loaded <- function(df) {
	plotter <- ggplot(df) + 
		geom_point(aes(x=af, y=tau, color=test_loaded),) +
             gtex_v8_figure_theme() + 
             labs(x="Simulated allele frequency", y = "residual precision", color="Test loaded strongest on strongest factor") + 
             theme(legend.position="bottom") + 
             theme(legend.text = element_text(size=8), legend.title = element_text(size=8))
  return(plotter)

}
make_af_vs_log_tau_scatter_plot_colored_by_whether_test_was_loaded <- function(df) {
	df$log_tau = log(df$tau)
	plotter <- ggplot(df) + 
		geom_point(aes(x=af, y=log_tau, color=test_loaded),) +
             gtex_v8_figure_theme() + 
             labs(x="Simulated allele frequency", y = "log(residual precision)", color="Test loaded strongest on strongest factor") + 
             theme(legend.position="bottom") + 
             theme(legend.text = element_text(size=8), legend.title = element_text(size=8))
  return(plotter)

}

make_tau_boxplot_stratefied_by_heteroskedasticity <- function(df) {
	boxplot <- ggplot(df, aes(x=version, y=tau)) + geom_boxplot() +
				gtex_v8_figure_theme() + 
	        	labs(x="", y = "Residual precision") +
	        	theme(legend.position="bottom")
	return(boxplot)
}

######################
# Command line args
######################
input_data_dir <- args[1]
visualization_dir <- args[2]


df <- load_in_data(input_data_dir)

output_file <- paste0(visualization_dir, "af_vs_tau_scatter_colored_by_test_loaded.pdf")
scatter <- make_af_vs_tau_scatter_plot_colored_by_whether_test_was_loaded(df[as.character(df$version) == "heteroskedasticity",] )
ggsave(scatter, file=output_file, width=7.2, height=5.5, units="in")

output_file <- paste0(visualization_dir, "af_vs_log_tau_scatter_colored_by_test_loaded.pdf")
scatter <- make_af_vs_log_tau_scatter_plot_colored_by_whether_test_was_loaded(df[as.character(df$version) == "heteroskedasticity",] )
ggsave(scatter, file=output_file, width=7.2, height=5.5, units="in")

output_file <- paste0(visualization_dir, "tau_boxplot_stratefied_by_heteroskedasticity.pdf")
boxplot <- make_tau_boxplot_stratefied_by_heteroskedasticity(df[df$af == .1,] )
ggsave(boxplot, file=output_file, width=7.2, height=5.5, units="in")


