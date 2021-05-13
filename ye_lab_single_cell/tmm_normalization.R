library(edgeR)









#########################
# Command line args
##########################
args <- commandArgs(TRUE)
input_file <- args[1]  


output_stem = strsplit(as.character(input_file), '[.]')[[1]][1]

# Load in input file
counts <- read.table(input_file, header=FALSE, sep="\t")
print(dim(counts))

# Create DGEList object
dgeFull <- DGEList(counts, remove.zeros = FALSE)
# Calculate TMM normalization factors
dgeFull <- calcNormFactors(dgeFull, method="TMMwzp")

# Generate CPMs
normCounts <- cpm(dgeFull)
# Generate log CPMs
pseudoNormCounts <- cpm(dgeFull, log=TRUE)

# Save cpm to output
output_stem = strsplit(as.character(input_file), '[.]')[[1]][1]
cpm_output_file <- paste0(output_stem, "_cpm.txt")
write.table(normCounts, cpm_output_file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

# Save log-cpm to output
log_cpm_output_file <- paste0(output_stem, "_log_cpm.txt")
write.table(pseudoNormCounts, log_cpm_output_file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
