library(scran)
library(scater)
library(Matrix)




##########################
# Command line args
##########################
args <- commandArgs(TRUE)
expr_file <- args[1] 
batch_file <- args[2] 

batch <- read.table(batch_file)

log_counts <- as.matrix(read.table(expr_file, header=FALSE, sep="\t"))
saveRDS(log_counts, "scran_input_log_counts.rds")
#log_counts <- readRDS("scran_input_log_counts.rds")

sce <- SingleCellExperiment(list(logcounts=log_counts))


design = model.matrix(~factor(batch$V1))
trend_fit <- trendVar(sce, design=design, use.spikes=FALSE)

dec <- decomposeVar(sce, trend_fit, design=design)


write.table(dec, "hvg.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")