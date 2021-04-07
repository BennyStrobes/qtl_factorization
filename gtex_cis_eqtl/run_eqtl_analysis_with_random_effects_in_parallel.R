args = commandArgs(trailingOnly=TRUE)
library(lme4)
library(lmtest)
library(ggplot2)
library(reshape2)
options(bitmapType='cairo')


run_eqtl_lmm <- function(expr, geno, covariates, groups) {
	fit_full <- lmer(expr ~ geno + covariates + (1|groups), REML=FALSE)
	fit_null <- lmer(expr ~ covariates + (1|groups), REML=FALSE)
	lrt <- anova(fit_null,fit_full)
	coefs <- data.frame(coef(summary(fit_full)))
	beta <- coefs[2,1]
	std_err <- coefs[2,2]
	pvalue <- lrt[[8]][2]
	return(list(eqtl_pvalue=pvalue, eqtl_beta=beta, eqtl_std_err=std_err))
}




#####################
# Command line args
#####################
expression_file = args[1]
genotype_file = args[2]
test_names_file = args[3]
covariate_file = args[4]
sample_overlap_file = args[5]
output_root = args[6]
job_number = as.numeric(args[7])
num_jobs = as.numeric(args[8])
total_lines = as.numeric(args[9])

print(total_lines)


# Determine number of lines each parrallelized job will complete
lines_per_job = ceiling(total_lines/num_jobs)
print(lines_per_job)
start_num = job_number*lines_per_job
end_num = (job_number + 1)*lines_per_job


covariates <- as.matrix(read.table(covariate_file, header=FALSE))
groups <- read.table(sample_overlap_file, header=FALSE)$V1


output_file <- paste0(output_root, "results.txt")
sink(output_file)

# Stream files
stop = FALSE
count = 0

f_test = file(test_names_file, "r")
f_expr = file(expression_file, "r")
f_geno = file(genotype_file, "r")


line_test = readLines(f_test, n=1)
line_test = readLines(f_test, n=1)
line_expr = readLines(f_expr, n=1)
line_geno = readLines(f_geno, n=1)

while(!stop) {
	if (count >= start_num & count < end_num) {
	# Unpack data
		expr = as.numeric(strsplit(line_expr,'\t')[[1]])
		geno = as.numeric(strsplit(line_geno,'\t')[[1]])
		# Run eqtl analysis
		line_info <- strsplit(line_test,'\t')[[1]]
		ensamble_id = line_info[1]
		rs_id = line_info[2]
		
		tryCatch(
		{
			lmm_results = run_eqtl_lmm(expr, geno, covariates, groups)
			new_line <- paste0(rs_id, "\t", ensamble_id ,"\t",lmm_results$eqtl_beta, "\t", lmm_results$eqtl_std_err,"\t", lmm_results$eqtl_pvalue,"\n")
        	cat(new_line)
        },
        error = function(e) {
        	new_line <- paste0(rs_id, "\t", ensamble_id ,"\t",0.0, "\t", 1.0,"\t", 1.0,"\n")
        }
        )
	}
	# Get info for next line
	count = count + 1
	line_test = readLines(f_test, n=1)
	line_expr = readLines(f_expr, n=1)
	line_geno = readLines(f_geno, n=1)
	if (length(line_test) == 0) {
		stop=TRUE
		close(f_test)
		close(f_expr)
		close(f_geno)
	}
}
sink()