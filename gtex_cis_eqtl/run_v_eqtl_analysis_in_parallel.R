args = commandArgs(trailingOnly=TRUE)
library(lme4)
library(lmtest)
library(ggplot2)
library(reshape2)
options(bitmapType='cairo')


run_v_eqtl_lm <- function(expr, geno, covariates) {
	fit_full <- lm(expr ~ geno + covariates)
	squared_geno <- geno*geno
	het_fit_full <- lm(abs(resid(fit_full)) ~ geno + squared_geno)
	het_fit_null <- lm(abs(resid(fit_full)) ~ 1)

    obj <- lrtest(het_fit_null, het_fit_full)
    pvalue <- obj[[5]][2]
    coefs <- data.frame(coef(summary(het_fit_full)))
    beta <- coefs[2,1]
    beta_squared <- coefs[3,1]
	return(list(eqtl_pvalue=pvalue, geno_coef=beta, geno_squared_coef=beta_squared))
}



#####################
# Command line args
#####################
expression_file = args[1]
genotype_file = args[2]
test_names_file = args[3]
covariate_file = args[4]
output_root = args[5]
job_number = as.numeric(args[6])
num_jobs = as.numeric(args[7])
total_lines = as.numeric(args[8])


# Determine number of lines each parrallelized job will complete
lines_per_job = ceiling(total_lines/num_jobs)
print(lines_per_job)
start_num = job_number*lines_per_job
end_num = (job_number + 1)*lines_per_job


covariates <- as.matrix(read.table(covariate_file, header=FALSE))


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
			lmm_results = run_v_eqtl_lm(expr, geno, covariates)
			new_line <- paste0(rs_id, "\t", ensamble_id ,"\t",lmm_results$geno_coef, "\t", lmm_results$geno_squared_coef,"\t", lmm_results$eqtl_pvalue,"\n")
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