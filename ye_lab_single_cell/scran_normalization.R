library(scran)


##########################
# Command line args
##########################
args <- commandArgs(TRUE)
input_file <- args[1]  

output_stem = strsplit(as.character(input_file), '[.]tx')[[1]][1]
print(output_stem)

# Load in input file
counts <- as.matrix(t(read.table(input_file, header=FALSE, sep="\t")))
#saveRDS(counts, "scran_input_raw_counts.rds")
#counts <- readRDS("scran_input_raw_counts.rds")



sce <- SingleCellExperiment(list(counts=counts))


clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clusters)
sce <- normalize(sce)


# Save log-cpm to output
log_cpm_output_file <- paste0(output_stem, "_log_norm_cpm.txt")
write.table(t(logcounts(sce)), log_cpm_output_file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

