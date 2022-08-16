args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
options(bitmapType = 'cairo', device = 'pdf')

figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


make_covariate_loading_correlation_heatmap <- function(covariates, loadings) {
	loadings <- as.matrix(loadings)
	covs <- covariates
    
    # Filter out covs that only contain 1 factor
    num_covs_init <- dim(covs)[2]
    cov_boolean <- c()
    for (num_cov in 1:num_covs_init) {
        cov_vec <- covs[, num_cov]
        if (is.factor(cov_vec) == TRUE) {
            if (nlevels(cov_vec) <= 1) {
                cov_boolean <- c(cov_boolean, FALSE)
            } else {
                cov_boolean <- c(cov_boolean, TRUE)
            }
        } else {
            cov_boolean <- c(cov_boolean, TRUE)
        }
    }
    covs <- covs[,cov_boolean]

    
    # Initialize PVE heatmap
    factor_colnames <- paste0("Factor", 1:(dim(loadings)[2]))
    factor_rownames <- colnames(covs)
    pve_map <- matrix(0, dim(covs)[2], dim(loadings)[2])
    colnames(pve_map) <- factor_colnames
    rownames(pve_map) <- colnames(covs)


    # Loop through each PC, COV Pair and take correlation
    num_pcs <- dim(loadings)[2]
    num_covs <- dim(covs)[2]
    for (num_pc in 1:num_pcs) {
        for (num_cov in 1:num_covs) {
            pc_vec <- loadings[,num_pc]
            cov_vec <- covs[,num_cov]
            lin_model <- lm(pc_vec ~ cov_vec)
            pve_map[num_cov, num_pc] <- summary(lin_model)$adj.r.squared
        }
    }
    pve_map[pve_map < 0] = 0.0

    ord <- hclust( dist(scale(pve_map), method = "euclidean"), method = "ward.D" )$order

    melted_mat <- melt(pve_map)
    colnames(melted_mat) <- c("Covariate", "Loading","PVE")

    melted_mat$Covariate = factor(melted_mat$Covariate, levels=rownames(pve_map)[ord])
    melted_mat$Loading = factor(melted_mat$Loading, levels=factor_colnames)
	 #  Use factors to represent covariate and pc name
   	# melted_mat$Covariate 
    # melted_mat$Covariate <- factor(melted_mat$Covariate, levels = rownames(pve_map)[ord])
    #melted_mat$PC <- substr(as.character(melted_mat$PC),3,5)
    #melted_mat$PC <- factor(melted_mat$PC, levels=paste0("", 1:(length(unique(melted_mat$PC)))))

    
    #levels(melted_mat$PC) = paste0("PC", 1:(length(levels(melted_mat$PC))))
    #  PLOT!
    heatmap <- ggplot(data=melted_mat, aes(x=Covariate, y=Loading)) + geom_tile(aes(fill=PVE)) + scale_fill_gradient2(midpoint=-.05, guide="colorbar")
    heatmap <- heatmap + labs(y="",fill="VE")
    heatmap <- heatmap + theme(text = element_text(size=7),axis.text=element_text(size=6), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8),  axis.text.x = element_text(angle = 90,hjust=1, vjust=.5)) 
    # Save File
    return(heatmap)
}

generate_per_tissue_covariate_df <- function(df) {
	df2 <- data.frame(df)
	unique_tissues <- as.character(unique(df2$tissue_type))
	num_tissues <- length(unique_tissues)
	for (tissue_num in 1:num_tissues) {
		tissue_name <- unique_tissues[tissue_num]
		tissue_indicator_vec <- 1.0*(as.character(df2$tissue_type) == tissue_name)
		df2[, tissue_name] <- tissue_indicator_vec
	}
	return(df2[,2:(dim(df2)[2])])
}


# Make plot showing variance explained of first n pcs
plot_pca_variance_explained <- function(sv, n) {

    variance_explained <- (sv$d^2/sum(sv$d^2))[1:n]

    # Merge global vectors into a data frame
    df <- data.frame(variance_explained = variance_explained, pc_num = 1:n)

    # PLOT AWAY
    line_plot <- ggplot(data=df, aes(x=pc_num, y=variance_explained)) +
                geom_line() +
                geom_point() +
                ylim(0,max(variance_explained) + .01) + 
                scale_x_continuous(breaks=seq(0, n, by=10)) +
                labs(x = "Number of PCs", y = "Variance Explained") + 
                theme(plot.title=element_text(size=8,face="plain"),text = element_text(size=8),axis.text=element_text(size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=7), legend.title = element_text(size=8), axis.text.x = element_text(vjust=.5)) 

    # SAVE PLOT
    return(line_plot)
}


tissues_file <- args[1]
data_dir <- args[2]
visualization_dir <- args[3]


# Load in Expression File
expr_file <- paste0(data_dir, "cross_tissue_tpm_standardized.txt")
expr <- read.table(expr_file, header=TRUE, sep="\t")
expr <- as.matrix(expr[,2:(dim(expr)[2])])
# Generate Expression PCs
svd1 <- svd(expr)
saveRDS(svd1, paste0(visualization_dir, "svd1.rds"))
# If already run (to save time)
#svd1 <- readRDS(paste0(visualization_dir, "svd1.rds"))


# Load in Expression File
#resid_expr_file <- paste0(data_dir, "cross_tissue_tpm_technical_covariate_residuals.txt")
#resid_expr <- read.table(resid_expr_file, header=TRUE, sep="\t")
#resid_expr <- as.matrix(resid_expr[,2:(dim(resid_expr)[2])])
# Generate Expression PCs
#svd2 <- svd(resid_expr)
#saveRDS(svd2, paste0(visualization_dir, "svd2.rds"))
# If already run (to save time)
#svd2 <- readRDS(paste0(visualization_dir, "svd2.rds"))

# Load in known covariates
individual_covariate_file <- paste0(data_dir, "sample_covariates.txt")
individual_covariates <- read.table(individual_covariate_file, header=TRUE, sep="\t")
individual_covariates <- individual_covariates[, 2:(dim(individual_covariates)[2])]
individual_covariates$race = factor(individual_covariates$race)

# Load in technical covariates
technical_covariate_file <- paste0(data_dir, "sample_technical_covariates.txt")
technical_covariates <- read.table(technical_covariate_file, header=TRUE, sep="\t")
technical_covariates <- technical_covariates[, 2:(dim(technical_covariates)[2])]

# Create individual covariate file with seperate column for each tissue
individual_per_tissue_covariates <- generate_per_tissue_covariate_df(individual_covariates)




######################################
# Make PCA variance explained line plot
#######################################
num_pcs <- 200
output_file <- paste0(visualization_dir, "expr_", num_pcs, "_pc_variance_explained_line_plot.pdf")
ve_lineplot <- plot_pca_variance_explained(svd1, num_pcs)
ggsave(ve_lineplot, file=output_file, width=7.2, height=6.0, units="in")

if (FALSE) {
num_pcs <- 200
output_file <- paste0(visualization_dir, "resid_expr_", num_pcs, "_pc_variance_explained_line_plot.pdf")
ve_lineplot <- plot_pca_variance_explained(svd2, num_pcs)
ggsave(ve_lineplot, file=output_file, width=7.2, height=6.0, units="in")
}
######################################
# Make correlation heatmap correlating covariates with loadings
#######################################
num_pcs <- 100
output_file <- paste0(visualization_dir, "expr_", num_pcs, "_pc_individual_covariates_correlation_heatmap.pdf")
heatmap <- make_covariate_loading_correlation_heatmap(individual_covariates, svd1$v[,1:num_pcs])
ggsave(heatmap, file=output_file, width=7.2, height=9.0, units="in")

if (FALSE) {
num_pcs <- 100
output_file <- paste0(visualization_dir, "resid_expr_", num_pcs, "_pc_individual_covariates_correlation_heatmap.pdf")
heatmap <- make_covariate_loading_correlation_heatmap(individual_covariates, svd2$v[,1:num_pcs])
ggsave(heatmap, file=output_file, width=7.2, height=9.0, units="in")
}

######################################
# Make correlation heatmap correlating covariates with loadings
#######################################
num_pcs <- 100
output_file <- paste0(visualization_dir, "expr_", num_pcs, "_pc_individual_per_tissue_covariates_correlation_heatmap.pdf")
#heatmap <- make_covariate_loading_correlation_heatmap(individual_per_tissue_covariates, svd1$v[,1:num_pcs])
#ggsave(heatmap, file=output_file, width=7.2, height=9.0, units="in")

num_pcs <- 100
output_file <- paste0(visualization_dir, "resid_expr_", num_pcs, "_pc_individual_per_tissue_covariates_correlation_heatmap.pdf")
#heatmap <- make_covariate_loading_correlation_heatmap(individual_per_tissue_covariates, svd2$v[,1:num_pcs])
#ggsave(heatmap, file=output_file, width=7.2, height=9.0, units="in")

######################################
# Make correlation heatmap correlating covariates with loadings
#######################################
num_pcs <- 100
output_file <- paste0(visualization_dir, "expr_", num_pcs, "_pc_sample_technical_covariates_correlation_heatmap.pdf")
heatmap <- make_covariate_loading_correlation_heatmap(technical_covariates, svd1$v[,1:num_pcs])
ggsave(heatmap, file=output_file, width=7.2, height=9.0, units="in")

if (FALSE) {
num_pcs <- 100
output_file <- paste0(visualization_dir, "resid_expr_", num_pcs, "_pc_sample_technical_covariates_correlation_heatmap.pdf")
heatmap <- make_covariate_loading_correlation_heatmap(technical_covariates, svd2$v[,1:num_pcs])
ggsave(heatmap, file=output_file, width=7.2, height=9.0, units="in")
}
