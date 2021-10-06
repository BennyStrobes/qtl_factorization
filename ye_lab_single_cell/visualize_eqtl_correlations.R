args = commandArgs(trailingOnly=TRUE)
library(reshape)
library(grid)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
options(bitmapType = 'cairo', device = 'pdf')


gtex_v8_figure_theme <- function() {
	return(theme(plot.title = element_text(face="plain",size=11), text = element_text(size=11),axis.text=element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.text = element_text(size=11), legend.title = element_text(size=11)))
}


symmetric_correlation_heatmap_general <- function(correlation_matrix) {
    nn <- dim(correlation_matrix)[1]
    vec <- c()
    for (i in 1:nn){
        for (j in 1:nn) {
            if (i != j) {
                vec <- c(vec,correlation_matrix[i,j])
            }
        }
    }
    maxy <- max(vec)
    for (i in 1:nn) {
        correlation_matrix[i,i] <- maxy
    }
    melted_corr <- melt(correlation_matrix)

    # Axis labels are factors
    melted_corr$X1 <- factor(melted_corr$X1, levels = rownames(correlation_matrix))
    melted_corr$X2 <- factor(melted_corr$X2, levels = colnames(correlation_matrix))

    #  PLOT!
    heatmap <- ggplot(data=melted_corr, aes(x=X1, y=X2)) + geom_tile(aes(fill=value)) 
    #heatmap <- heatmap + scale_fill_distiller()
    #heatmap <- heatmap + scale_fill_brewer(values = brewer.pal(3,"RdPu"))
    heatmap <- heatmap + scale_fill_distiller(palette = "RdPu", direction=1)
    heatmap <- heatmap + theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90))
    heatmap <- heatmap + gtex_v8_figure_theme()
    heatmap <- heatmap + labs(fill= "Spearman Rho",x="", y="")

    return(heatmap)
}







eqtl_correlation_dir <- args[1]

# Correlation of SURGE context specific loadings
eqtl_loading_file <- paste0(eqtl_correlation_dir, "loading_correlation.txt")
loading_correlation_mat <- as.matrix(read.table(eqtl_loading_file, header=FALSE, sep="\t"))
rownames(loading_correlation_mat) = paste0("surge_context_", 1:10)
colnames(loading_correlation_mat) = paste0("surge_context_", 1:10)

output_file <- paste0(eqtl_correlation_dir, "surge_loading_correlation_heatmap.pdf")
heatmap <- symmetric_correlation_heatmap_general(loading_correlation_mat)
ggsave(heatmap, file=output_file, width=7.2, height=6.0, units="in")


# Correlation of SURGE context specific eqtl summary stats
eqtl_correlation_file <- paste0(eqtl_correlation_dir, "eqtl_correlation.txt")
eqtl_correlation_mat <- as.matrix(read.table(eqtl_correlation_file, header=FALSE, sep="\t"))
rownames(eqtl_correlation_mat) = paste0("surge_context_", 1:10, "_eqtl")
colnames(eqtl_correlation_mat) = paste0("surge_context_", 1:10, "_eqtl")

output_file <- paste0(eqtl_correlation_dir, "eqtl_correlation_heatmap.pdf")
heatmap <- symmetric_correlation_heatmap_general(eqtl_correlation_mat)
ggsave(heatmap, file=output_file, width=7.2, height=6.0, units="in")

# Correlation of SURGE context specific eqtl summary stats + standard eqtls
eqtl_correlation_file <- paste0(eqtl_correlation_dir, "all_eqtl_correlation.txt")
eqtl_correlation_mat <- as.matrix(read.table(eqtl_correlation_file, header=FALSE, sep="\t"))
rownames(eqtl_correlation_mat) = c(paste0("surge_context_", 1:10, "_eqtl"), "standard_eqtl")
colnames(eqtl_correlation_mat) = c(paste0("surge_context_", 1:10, "_eqtl"), "standard_eqtl")

output_file <- paste0(eqtl_correlation_dir, "all_eqtl_correlation_heatmap.pdf")
heatmap <- symmetric_correlation_heatmap_general(eqtl_correlation_mat)
ggsave(heatmap, file=output_file, width=7.2, height=6.0, units="in")
