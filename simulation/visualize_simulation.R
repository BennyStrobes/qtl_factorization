args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(cowplot)
library(umap)
library(ggplot2)
library(RColorBrewer)
library(lme4)
library(sigmoid)
library(ggbeeswarm)
options(bitmapType = 'cairo', device = 'pdf')



figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


make_umap_loading_scatter_plot_colored_by_real_valued_variable <- function(covariates, umap_loadings, covariate_name,point_size=1) {
	df <- data.frame(loading_1=umap_loadings[,1], loading_2=umap_loadings[,2], covariate=covariates)
	plotter <- ggplot(df) + 
	           geom_point(aes(x=loading_1, y=loading_2, color=covariate), size=point_size) +
	           figure_theme() + 
	           labs(x="", y = "", color=covariate_name) + 
               scale_colour_gradient2() +
	           scale_color_gradient(low="pink",high="blue") +
	           theme(legend.text = element_text(size=8), legend.title = element_text(size=8)) +
	           theme_void() +
	           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank()) 
	return(plotter)
}

make_simulated_interaction_eqtl_plot <- function(Y, G, context, context_legend, point_size=.65) {
	G_str = rep("AA", length(G))
	G_str[G==1] = "AT"
	G_str[G==2] = "TT"
	df <- data.frame(expression=Y, genotype=factor(G_str), context=context)
	plotter <- ggplot(df) + 
	           geom_point(aes(x=context, y=expression, color=genotype), size=point_size) +
	           figure_theme() + 
	           labs(x="\n\nCellular context", y = "Gene Expression", color="") + 
	           theme(legend.position="right") + 
	           theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
	           theme(legend.text = element_text(size=8), legend.title = element_text(size=8))
	return(plotter)
}

make_ve_swarm_plot <- function(df) {
	df$t.statistic = factor(df$t.statistic)
	df$sample_size = factor(df$sample_size)
	plotter <- ggplot(df) +  geom_beeswarm(aes(x=sample_size, y=average_r_squared, color=t.statistic),size=.4) +
	figure_theme() +
	theme(legend.text = element_text(size=8), legend.title = element_text(size=11)) +
	labs(x="Sample size", y="Variance explained\nof simulated context", color="Interaction\nvariance")
	return(plotter)
}



make_num_contexts_boxplot <- function(df) {
	#df$number_simulated_components = factor(df$number_simulated_components)
	#df$number_learned_components = factor(df$number_learned_components)

	plotter <- ggplot(df) +  geom_beeswarm(aes(x=number_simulated_components, y=number_learned_components), color='skyblue', size=.4) +
	figure_theme() +
	scale_y_continuous(breaks=seq(0,21,2)) +
	scale_x_continuous(breaks=seq(0,21,2)) +
	geom_abline(size=.1,color='grey')+
	theme(legend.text = element_text(size=10), legend.title = element_text(size=11)) +
	labs(x="Number of simulated contexts", y="Number of learned contexts")
	return(plotter)
}

make_run_time_eqtl_lmm_analysis_boxplot <- function(eqtl_results_dir) {

	eqtl_ss_vec <- c()
	run_time_vec <- c()
	sample_vec <- c()
	for (itera in 1:1) {
		tmp_df <- read.table(paste0(eqtl_results_dir, "eqtl_lmm_run_time_", itera, ".txt"), header=TRUE)
		#print(tmp_df)
		tmp_df = tmp_df[2:dim(tmp_df)[1],]

		eqtl_ss_vec <- c(eqtl_ss_vec, tmp_df$sample_size)
		run_time_vec <- c(run_time_vec, tmp_df$run_time)
	}

	df <- data.frame(eqtl_ss=factor(eqtl_ss_vec), run_time=run_time_vec)
	pp <- ggplot(df, aes(x=eqtl_ss, y=run_time)) + 
 		 geom_beeswarm(color='blue', size=.5, cex=.54) + 
 		 figure_theme() +
 		 labs(x="eQTL sample size", y="Run time (seconds)", title="LME4 Interaction eQTL Runtime")

 	return(pp)
}

make_run_time_analysis_boxplot <- function(eqtl_results_dir) {
	eqtl_ss_vec <- c()
	run_time_vec <- c()
	sample_vec <- c()
	for (itera in 1:10) {
		tmp_df <- read.table(paste0(eqtl_results_dir, "run_time_analysis_", itera, ".txt"), header=TRUE)
		tmp_df = tmp_df[2:dim(tmp_df)[1],]

		eqtl_ss_vec <- c(eqtl_ss_vec, tmp_df$eqtl_sample_size)
		run_time_vec <- c(run_time_vec, tmp_df$run_time/(60.0*60.0))
		sample_vec <- c(sample_vec, tmp_df$sample)
	}

	df <- data.frame(eqtl_ss=factor(eqtl_ss_vec), run_time=run_time_vec)
	pp <- ggplot(df, aes(x=eqtl_ss, y=run_time)) + 
 		 geom_beeswarm(color='blue', size=1.0) + 
 		 figure_theme() +
 		 labs(x="eQTL sample size", y="Run time (hours)", title="SURGE optimization runtime")

 	return(pp)	
}


eqtl_results_dir=args[1]
viz_dir=args[2]


################################################
# Make run time analysis boxplot
################################################
output_file <- paste0(viz_dir, "surge_runtime_boxplot.pdf")
boxplot <- make_run_time_analysis_boxplot(eqtl_results_dir)
ggsave(boxplot, file=output_file, width=7.2, height=4.0, units="in")

output_file <- paste0(viz_dir, "eqtl_lmm_runtime_boxplot.pdf")
boxplot2 <- make_run_time_eqtl_lmm_analysis_boxplot(eqtl_results_dir)
ggsave(boxplot2, file=output_file, width=7.2, height=4.0, units="in")

# make joint runtime plot
joint_plot <- plot_grid(boxplot, boxplot2, ncol=1)
output_file <- paste0(viz_dir, "joint_runtime.pdf")
ggsave(joint_plot, file=output_file, width=7.2, height=5.5, units="in")


if (FALSE) {
################################################
# Assess ability to identify correct number of simulated factors
################################################
num_contexts_1_3_results_file <- paste0(eqtl_results_dir, "number_of_learned_components_recovery_results_250_0.1_0.3.txt")
num_contexts_1_3_results <- read.table(num_contexts_1_3_results_file, header=TRUE)

num_contexts_25_3_results_file <- paste0(eqtl_results_dir, "number_of_learned_components_recovery_results_250_0.25_0.3.txt")
num_contexts_25_3_results <- read.table(num_contexts_25_3_results_file, header=TRUE)

num_contexts_5_3_results_file <- paste0(eqtl_results_dir, "number_of_learned_components_recovery_results_250_0.5_0.3.txt")
num_contexts_5_3_results <- read.table(num_contexts_5_3_results_file, header=TRUE)

num_contexts_1_1_results_file <- paste0(eqtl_results_dir, "number_of_learned_components_recovery_results_250_0.1_0.1.txt")
num_contexts_1_1_results <- read.table(num_contexts_1_1_results_file, header=TRUE)

num_contexts_25_1_results_file <- paste0(eqtl_results_dir, "number_of_learned_components_recovery_results_250_0.25_0.1.txt")
num_contexts_25_1_results <- read.table(num_contexts_25_1_results_file, header=TRUE)

num_contexts_5_1_results_file <- paste0(eqtl_results_dir, "number_of_learned_components_recovery_results_250_0.5_0.1.txt")
num_contexts_5_1_results <- read.table(num_contexts_5_1_results_file, header=TRUE)

output_file <- paste0(viz_dir, "simulation_num_contexts_boxplot_1_1.pdf")
num_contexts_boxplot_1_1 <- make_num_contexts_boxplot(num_contexts_1_1_results)
ggsave(num_contexts_boxplot_1_1, file=output_file, width=7.2, height=4.0, units="in")

output_file <- paste0(viz_dir, "simulation_num_contexts_boxplot_25_1.pdf")
num_contexts_boxplot_25_1 <- make_num_contexts_boxplot(num_contexts_25_1_results)
ggsave(num_contexts_boxplot_25_1, file=output_file, width=7.2, height=4.0, units="in")

output_file <- paste0(viz_dir, "simulation_num_contexts_boxplot_5_1.pdf")
num_contexts_boxplot_5_1 <- make_num_contexts_boxplot(num_contexts_5_1_results)
ggsave(num_contexts_boxplot_5_1, file=output_file, width=7.2, height=4.0, units="in")


output_file <- paste0(viz_dir, "simulation_num_contexts_boxplot_1_3.pdf")
num_contexts_boxplot_1_3 <- make_num_contexts_boxplot(num_contexts_1_3_results)
ggsave(num_contexts_boxplot_1_3, file=output_file, width=7.2, height=4.0, units="in")

output_file <- paste0(viz_dir, "simulation_num_contexts_boxplot_25_3.pdf")
num_contexts_boxplot_25_3 <- make_num_contexts_boxplot(num_contexts_25_3_results)
ggsave(num_contexts_boxplot_25_3, file=output_file, width=7.2, height=4.0, units="in")

output_file <- paste0(viz_dir, "simulation_num_contexts_boxplot_5_3.pdf")
num_contexts_boxplot_5_3 <- make_num_contexts_boxplot(num_contexts_5_3_results)
ggsave(num_contexts_boxplot_5_3, file=output_file, width=7.2, height=4.0, units="in")


output_file <- paste0(viz_dir, "supplement_simulation_num_contexts_boxplot_joint.pdf")
joint_num_contexts <- plot_grid(num_contexts_boxplot_25_1 +labs(title="Interaction term variance: .25\nFraction of non-zero tests: .1") , num_contexts_boxplot_25_3+labs(title="Interaction term variance: .25\nFraction of non-zero tests: .3"), num_contexts_boxplot_5_1+labs(title="Interaction term variance: .5\nFraction of non-zero tests: .1"), num_contexts_boxplot_5_3+labs(title="Interaction term variance: .5\nFraction of non-zero tests: .3"), labels=c("A","B","C","D"), ncol=2)
ggsave(joint_num_contexts, file=output_file, width=7.2, height=6.5, units="in")



################################################
# Assess Fraction of variance explained by simulated latent factors
################################################
variance_explained_1_results_file <- paste0(eqtl_results_dir, "variance_explained_power_analysis_missingness_0.1_results.txt")
variance_explained_1_df <- read.table(variance_explained_1_results_file, header=TRUE)
output_file <- paste0(viz_dir, "simulation_variance_explained_1.pdf")
ve_swarm_plot_1 <- make_ve_swarm_plot(variance_explained_1_df)
ggsave(ve_swarm_plot_1, file=output_file, width=7.2, height=4.0, units="in")

variance_explained_3_results_file <- paste0(eqtl_results_dir, "variance_explained_power_analysis_missingness_0.3_results.txt")
variance_explained_3_df <- read.table(variance_explained_3_results_file, header=TRUE)
output_file <- paste0(viz_dir, "simulation_variance_explained_3.pdf")
ve_swarm_plot_3 <- make_ve_swarm_plot(variance_explained_3_df)
ggsave(ve_swarm_plot_3, file=output_file, width=7.2, height=4.0, units="in")

variance_explained_5_results_file <- paste0(eqtl_results_dir, "variance_explained_power_analysis_missingness_0.5_results.txt")
variance_explained_5_df <- read.table(variance_explained_5_results_file, header=TRUE)
output_file <- paste0(viz_dir, "simulation_variance_explained_5.pdf")
ve_swarm_plot_5 <- make_ve_swarm_plot(variance_explained_5_df)
ggsave(ve_swarm_plot_5, file=output_file, width=7.2, height=4.0, units="in")


#fig_1C <- make_ve_swarm_plot(variance_explained_3_df[variance_explained_3_df$sample_size!=1000.0,])
fig_1C <- make_ve_swarm_plot(variance_explained_3_df)

# Make joint plot
output_file <- paste0(viz_dir, "supplement_simulation_variance_explained.pdf")
joint_ve_plot = plot_grid(ve_swarm_plot_1 + labs(title="Fraction of variant-gene pairs with non-zero genotype interaction: .1"), ve_swarm_plot_3+ labs(title="Fraction of variant-gene pairs with non-zero genotype interaction: .3"), ve_swarm_plot_5+ labs(title="Fraction of variant-gene pairs with non-zero genotype interaction: .5"), labels=c("A","B","C"), ncol = 1)
ggsave(joint_ve_plot, file=output_file, width=7.2, height=7.0, units="in")



output_file <- paste0(viz_dir, "Fig_1CD.pdf")
fig1CD_tmp = plot_grid(fig_1C + theme(legend.position="bottom"), num_contexts_boxplot_25_3, nrow = 1, labels = c('C', 'D'))
ggsave(fig1CD_tmp, file=output_file, width=7.2, height=3.0, units="in")





################################################
# Generate schematic UMAP plot from simulated data
################################################
set.seed(1) 
N=3000
T=500
K=1
K2=6
U = matrix( rnorm(N*K,mean=0,sd=1), N, 1)
V = matrix( rnorm(K*T,mean=0,sd=1), 1, T)
U2 = matrix( rnorm(N*K2,mean=0,sd=1), N, K2)
V2 = matrix( rnorm(K2*T,mean=0,sd=1), K2, T)

U3 = matrix(sample(c(0,1), replace=TRUE, prob=c(.7,.3), size=N), N,1)
V3 = matrix( rnorm(1*T,mean=0,sd=1), 1, T)


mean_expr <- U%*%V + .4*(U2%*%V2) + U3%*%V3*.6
expr_noise = matrix( rnorm(N*T,mean=0,sd=.01), N, T)
expr = mean_expr + expr_noise

expr_umap = umap(expr)$layout


output_file <- paste0(viz_dir, "simulated_expression_umap.pdf")
sim_umap_plot = make_umap_loading_scatter_plot_colored_by_real_valued_variable(U[,1], expr_umap, "Cellular context", point_size=.1)
ggsave(sim_umap_plot, file=output_file, width=7.2, height=6.0, units="in")

context_legend <- get_legend(sim_umap_plot + labs(color="") + theme(legend.position="bottom") + theme(legend.key.height= unit(.26, 'cm'),legend.key.width= unit(.96, 'cm')))



################################################
# Generate schematic GXE plot from simulated data
################################################
G = sample(c(0,1), replace=TRUE, prob=c(.7,.3), size=N) + sample(c(0,1), replace=TRUE, prob=c(.7,.3), size=N)
mean_Y = G*(U[,1] - min(U[,1]) +.8)*.3
noise_Y = rnorm(N,mean=0,sd=.3)
Y = mean_Y + noise_Y

output_file <- paste0(viz_dir, "simulated_interaction_eqtl.pdf")
sim_interaction_eqtl_plot = make_simulated_interaction_eqtl_plot(Y, G, U[,1], context_legend)
ggsave(sim_interaction_eqtl_plot, file=output_file, width=7.2, height=4.0, units="in")



output_file <- paste0(viz_dir, "Fig_1A.pdf")
fig1a_tmp = plot_grid(sim_umap_plot, NULL,sim_interaction_eqtl_plot, nrow = 1, rel_widths=c(.8,.15,1), align='v')
fig1a <- ggdraw() + 
				draw_plot(fig1a_tmp,0,0,1,1) +
				draw_plot(context_legend,.185,-.29,1,1) 
ggsave(plot_grid(fig1a, labels = c("A")), file=output_file, width=7.2, height=2.2, units="in")

}




