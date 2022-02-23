args = commandArgs(trailingOnly=TRUE)
library(lme4)
library(lmtest)
library(ggplot2)
library(reshape2)
options(bitmapType='cairo')


run_eqtl_lm <- function(expr, geno, covariates) {
	fit_full <- lm(expr ~ geno + covariates)
	coefs <- data.frame(coef(summary(fit_full)))
	beta <- coefs[2,1]
	std_err <- coefs[2,2]
	pvalue <- coefs[2,4]
	return(list(eqtl_pvalue=pvalue, eqtl_beta=beta, eqtl_std_err=std_err))
}
run_eqtl_lmm <- function(expr, geno, covariates, groups) {
	fit_full <- lmer(expr ~ geno + covariates + (1 | groups), REML=FALSE)


	# extract coefficients
	coefs <- data.frame(coef(summary(fit_full)))
	effect_size <- coefs[2,1]
	std_err <- coefs[2,2]
	t_value <- coefs[2,3]
	aggregate_pvalue <- 2 * (1 - pnorm(abs(t_value)))



	return(list(eqtl_pvalue=aggregate_pvalue, eqtl_beta=effect_size, eqtl_std_err=std_err))
}




#####################
# Command line args
#####################
expression_file = args[1]
genotype_file = args[2]
test_names_file = args[3]
covariate_file = args[4]
qtl_sample_overlap_file = args[5]
output_root = args[6]
version = args[7]
job_number = as.numeric(args[8])
num_jobs = as.numeric(args[9])
total_lines = as.numeric(args[10]) -1
print(version)


# Determine number of lines each parrallelized job will complete
lines_per_job = ceiling(total_lines/num_jobs)
print(lines_per_job)
start_num = job_number*lines_per_job
end_num = (job_number + 1)*lines_per_job


covariates <- as.matrix(read.table(covariate_file, header=FALSE))
groups <- read.table(qtl_sample_overlap_file, header=FALSE)$V1 +1



output_file <- paste0(output_root, version, "_results.txt")
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

count = 0
while(!stop) {
	if (count >= start_num & count < end_num) {
	# Unpack data
		expr = as.numeric(strsplit(line_expr,'\t')[[1]])
		geno = as.numeric(strsplit(line_geno,'\t')[[1]])
		norm_geno = (geno - mean(geno))/(sd(geno))
		# Run eqtl analysis
		line_info <- strsplit(line_test,'\t')[[1]]
		ensamble_id = line_info[1]
		rs_id = line_info[2]

		tryCatch(
		{
			if (version == "lm") {
				lmm_results = run_eqtl_lm(expr, norm_geno, covariates)
			} else if (version == "lmm") {
				lmm_results = run_eqtl_lmm(expr, norm_geno, covariates, groups)
			}
			new_line <- paste0(rs_id, "\t", ensamble_id ,"\t",lmm_results$eqtl_beta, "\t", lmm_results$eqtl_std_err,"\t", lmm_results$eqtl_pvalue,"\n")
        	cat(new_line)
        },
        error = function(e) {
        	new_line <- paste0(rs_id, "\t", ensamble_id ,"\t",0.0, "\t", 1.0,"\t", 1.0,"\n")
        	cat(new_line)
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