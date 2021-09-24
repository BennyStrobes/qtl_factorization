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

run_lf_interaction_eqtl_lm <- function(expr, geno, covariates, lfs) {
	fit_full <- lm(expr ~ geno + covariates + lfs:geno)
	fit_null <- lm(expr ~ geno + covariates)
	obj <- lrtest(fit_null, fit_full)
	#lrt2 <- anova(fit_null, fit_full)
	aggregate_pvalue <- obj[[5]][2]
	num_cov = dim(covariates)[2]
	coefficient_pvalues = data.frame(coef(summary(fit_full)))[,4]
	lf_interaction_coefficient_pvalues = coefficient_pvalues[(3+num_cov):length(coefficient_pvalues)]
	return(list(eqtl_pvalue=aggregate_pvalue, coefficient_pvalues=lf_interaction_coefficient_pvalues))
}

run_lf_interaction_eqtl_lm_perm <- function(expr, geno, covariates, lfs, groups, individual_groups) {
	individual_genotype <- geno[individual_groups]
	if (all.equal(geno, individual_genotype[groups]) == FALSE) {
		print("EROROROR")
	}
	perm_geno <- (sample(individual_genotype))[groups]

	fit_full <- lm(expr ~ geno + covariates + lfs:perm_geno)
	fit_null <- lm(expr ~ geno + covariates)
	obj <- lrtest(fit_null, fit_full)
	#lrt2 <- anova(fit_null, fit_full)
	aggregate_pvalue <- obj[[5]][2]
	num_cov = dim(covariates)[2]
	coefficient_pvalues = data.frame(coef(summary(fit_full)))[,4]
	lf_interaction_coefficient_pvalues = coefficient_pvalues[(3+num_cov):length(coefficient_pvalues)]
	return(list(eqtl_pvalue=aggregate_pvalue, coefficient_pvalues=lf_interaction_coefficient_pvalues))
}

run_lf_interaction_eqtl_lmm_perm <- function(expr, geno, covariates, lfs, groups, individual_groups) {
	individual_genotype <- geno[individual_groups]
	if (all.equal(geno, individual_genotype[groups]) == FALSE) {
		print("EROROROR")
	}
	perm_geno <- (sample(individual_genotype))[groups]
	#fit_full <- lmer(expr ~ geno + covariates + lfs:geno + (1 | groups) + (0+lfs|groups), REML=FALSE)
	#fit_null <- lmer(expr ~ geno + covariates + (1 | groups) + (0+lfs|groups), REML=FALSE)

	fit_full <- lmer(expr ~ geno + covariates + lfs:perm_geno + (1 | groups), REML=FALSE)
	fit_null <- lmer(expr ~ geno + covariates + (1 | groups), REML=FALSE)

	lrt <- anova(fit_null,fit_full)

	aggregate_pvalue <- lrt[[8]][2]
	#obj <- lrtest(fit_null, fit_full)
	#aggregate_pvalue <- obj[[5]][2]
	#num_cov = dim(covariates)[2]
	#coefficient_pvalues = data.frame(coef(summary(fit_full)))[,4]
	#lf_interaction_coefficient_pvalues = coefficient_pvalues[(3+num_cov):length(coefficient_pvalues)]
	return(list(eqtl_pvalue=aggregate_pvalue, coefficient_pvalues=c(1.0, 1.0)))
}

run_lf_interaction_eqtl_lmm <- function(expr, geno, covariates, lfs, groups) {
	fit_full <- lmer(expr ~ geno + covariates + lfs + lfs:geno + (1 | groups), REML=FALSE)
	#fit_null <- lmer(expr ~ geno + covariates + lfs + (1 | groups), REML=FALSE)

	#lrt <- anova(fit_null,fit_full)

	#aggregate_pvalue <- lrt[[8]][2]

	# extract coefficients
	coefs <- data.frame(coef(summary(fit_full)))
	# use normal distribution to approximate p-value
	coefficient_pvalues <- 2 * (1 - pnorm(abs(coefs$t.value)))

	num_cov = dim(covariates)[2]
	num_lf = dim(lfs)[2]
	lf_interaction_coefficient_pvalues = coefficient_pvalues[(3+num_cov + num_lf):length(coefficient_pvalues)]
	betas <- coefs[(3+num_cov + num_lf):length(coefficient_pvalues),1]
	std_err <- coefs[(3+num_cov + num_lf):length(coefficient_pvalues),2]

	return(list(pvalue=lf_interaction_coefficient_pvalues, beta=betas, std_err=std_err))
}


pass_genotype_filter <- function(geno, thresh) {
	rounded_geno <- round(geno)
	booly = TRUE
	num_samples = length(rounded_geno)

	num_zero = sum(rounded_geno == 0)
	num_one = sum(rounded_geno == 1)
	num_two = sum(rounded_geno == 2)

	fraction_zero = num_zero/num_samples
	fraction_one = num_one/num_samples
	fraction_two = num_two/num_samples

	if (fraction_zero < thresh) {
		booly = FALSE
	}
	if (fraction_one < thresh) {
		booly = FALSE
	}
	if (fraction_two < thresh) {
		booly = FALSE
	}
	return(booly)
}

#####################
# Command line args
#####################
expression_file = args[1]
genotype_file = args[2]
test_names_file = args[3]
covariate_file = args[4]
interaction_factor_file = args[5]
sample_overlap_file = args[6]
output_root = args[7]
job_number = as.numeric(args[8])
num_jobs = as.numeric(args[9])
total_lines = as.numeric(args[10]) - 1

print(total_lines)


# Determine number of lines each parrallelized job will complete
lines_per_job = ceiling(total_lines/num_jobs)
print(lines_per_job)
start_num = job_number*lines_per_job
end_num = (job_number + 1)*lines_per_job


covariates <- as.matrix(read.table(covariate_file, header=FALSE))
lfs <- as.matrix(read.table(interaction_factor_file, header=FALSE))
groups <- read.table(sample_overlap_file, header=FALSE)$V1 +1


print("Data loaded in.. starting")

num_lfs <- dim(lfs)[2]
print(num_lfs)

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
			lmm_results = run_lf_interaction_eqtl_lmm(expr, geno, covariates, lfs, groups)

			new_line <- paste0(rs_id, "\t", ensamble_id ,"\t",lmm_results$beta, "\t", lmm_results$std_err, "\t", paste0(lmm_results$pvalue, collapse=","), "\n")
        		
        	cat(new_line)
        },
        error = function(e) {
        	new_line <- paste0(rs_id, "\t", ensamble_id, "\t", 0.0 ,"\t", 1.0, "\t", paste0(rep(1.0, num_lfs), collapse=","), "\n")
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