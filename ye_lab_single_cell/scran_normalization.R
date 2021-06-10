library(scran)
library(scater)
library(Matrix)
##########################
# Command line args
##########################
args <- commandArgs(TRUE)
input_file <- args[1]  

output_stem = strsplit(as.character(input_file), '[.]tx')[[1]][1]
print(output_stem)

# Load in input file
#counts <- as.matrix(read.table(input_file, header=FALSE, sep="\t"))
counts <- readSparseCounts(input_file, sep="\t", row.names=FALSE, col.names=FALSE)
print(dim(counts))
print('loaded')
#saveRDS(counts, "scran_input_raw_counts.rds")
#counts <- readRDS("scran_input_raw_counts.rds")



sce <- SingleCellExperiment(list(counts=counts))

print("A")
clusters <- quickCluster(sce)
print("B")
sce <- computeSumFactors(sce, cluster=clusters)
print("C")
sce <- normalize(sce)
print("D")


# Save log-cpm to output
log_cpm_output_file <- paste0(output_stem, "_log_norm_cpm_matrix_market_format.txt")
#write.table(as.matrix(t(logcounts(sce))), log_cpm_output_file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
writeMM(logcounts(sce), log_cpm_output_file)
