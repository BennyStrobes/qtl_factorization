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

run_lf_interaction_eqtl_lm <- function(expr, geno, perm_geno, covariates, lfs, num_lf) {
	fit_full <- lm(expr ~ geno + covariates + lfs + lfs:perm_geno)

	num_cov = dim(covariates)[2]

	summary_df <- data.frame(coef(summary(fit_full)))

	coefficient_betas = summary_df[,1]
	coefficient_std_err = summary_df[,2]
	coefficient_pvalues = summary_df[,4]

	lf_interaction_coefficient = coefficient_betas[(3+num_cov+num_lf):length(coefficient_betas)]
	lf_interaction_std_err = coefficient_std_err[(3+num_cov+num_lf):length(coefficient_std_err)]
	lf_interaction_pvalue = coefficient_pvalues[(3+num_cov+num_lf):length(coefficient_pvalues)]

	marginal_beta = coefficient_betas[2]
	marginal_std_err = coefficient_std_err[2]
	marginal_pvalue = coefficient_pvalues[2]

	return(list(pvalue=lf_interaction_pvalue, beta=lf_interaction_coefficient, std_err=lf_interaction_std_err, marginal_beta=marginal_beta, marginal_std_err=marginal_std_err, marginal_pvalue=marginal_pvalue))

}

get_huber_se <- function(H, hh_inv, fit_full, Y, N, little_N) {
	
	#little_N = 1000
	#little_N = 20000
	random_subsample = sample(nrow(H),size=little_N,replace=FALSE)
	
	squared_resid = (Y - predict(fit_full))^2
	Sigma = diag(squared_resid[random_subsample])


	H_sub = H[random_subsample,]

	inner = (1/little_N)*t(H_sub)%*%Sigma%*%H_sub

	full = hh_inv%*%inner%*%hh_inv

	huber_se = sqrt((1/N)*diag(full))
	return(huber_se)
}


## Extract the conditional variance-covariance matrix of the fixed-effects
## parameters
vcov.merMod <- function(object, correlation = TRUE, sigm = sigma(object),
                        use.hessian = NULL, ...)
{

    hess.avail <-
         ## (1) numerical Hessian computed?
        (!is.null(h <- object@optinfo$derivs$Hessian) &&
         ## (2) does Hessian include fixed-effect parameters?
         nrow(h) > (ntheta <- length(getME(object,"theta"))))
    if (is.null(use.hessian)) use.hessian <- hess.avail
    if (use.hessian && !hess.avail) stop(shQuote("use.hessian"),
                                         "=TRUE specified, ",
                                         "but Hessian is unavailable")
    calc.vcov.hess <- function(h) {
        ## invert 2*Hessian, catching errors and forcing symmetric result
        ## ~= forceSymmetric(solve(h/2)[i,i]) : solve(h/2) = 2*solve(h)
        h <- tryCatch(solve(h),
                      error=function(e) matrix(NA,nrow=nrow(h),ncol=ncol(h)))
        i <- -seq_len(ntheta)  ## drop var-cov parameters
        h <- h[i,i]
        forceSymmetric(h + t(h))
    }

    ## alternately: calculate var-cov from implicit (RX) information
    ## provided by fit (always the case for lmerMods)
    V <- sigm^2 * object@pp$unsc()

    if (hess.avail) {
        V.hess <- calc.vcov.hess(h)
        bad.V.hess <- any(is.na(V.hess))
        if (!bad.V.hess) {
            ## another 'bad var-cov' check: positive definite?
            e.hess <- eigen(V.hess,symmetric = TRUE,only.values = TRUE)$values
            if (min(e.hess) <= 0) bad.V.hess <- TRUE
        }
    }
    if (!use.hessian && hess.avail) {
        ## if hessian is available, go ahead and check
        ## for similarity with the RX-based estimate
        var.hess.tol <- 1e-4 # FIXME: should var.hess.tol be user controlled?
        if (!bad.V.hess && any(abs(V-V.hess) > var.hess.tol * V.hess))
            warning("variance-covariance matrix computed ",
                    "from finite-difference Hessian\nand ",
                    "from RX differ by >",var.hess.tol,": ",
                    "consider ",shQuote("use.hessian=TRUE"))
    }
    if (use.hessian) {
        if (!bad.V.hess) {
            V <- V.hess
        } else {
            warning("variance-covariance matrix computed ",
                    "from finite-difference Hessian is\n",
                    "not positive definite or contains NA values: falling back to ",
                    "var-cov estimated from RX")
        }
    }

    ## FIXME: try to catch non-PD matrices
    rr <- tryCatch(as(V, "dpoMatrix"), error = function(e)e)
    if (inherits(rr, "error")) {
        warning(gettextf("Computed variance-covariance matrix problem: %s;\nreturning NA matrix",
                         rr$message), domain = NA)
        rr <- matrix(NA,nrow(V),ncol(V))
    }

    nmsX <- colnames(object@pp$X)
    dimnames(rr) <- list(nmsX,nmsX)

    if(correlation)
        rr@factors$correlation <-
            if(!is.na(sigm)) as(rr, "corMatrix") else rr # (is NA anyway)
    rr
}


run_lf_interaction_eqtl_lmm <- function(expr, geno, perm_geno, covariates, lfs, groups, num_lf) {
	fit_full <- lmer(expr ~ geno + covariates + lfs + lfs:perm_geno + (1 | groups), REML=FALSE)
	#print(sqrt(diag(summary(fit_full)$vcov)))
	#print(((summary(fit_full)$sigma)))
	#print(summary(fit_full))
	N = length(expr)
	#hh_inv = as.matrix(vcov.merMod(fit_full))

	H=cbind(1, geno, covariates, lfs, lfs*perm_geno)

	#hh_inv = solve((1/N)*as.matrix(t(H)%*%H))
	#hh_inv = (N)*solve(as.matrix(t(H)%*%H))
	sigma_squared = (summary(fit_full)$sigma)^2
	vcov = vcov.merMod(fit_full)
	#hh_inv3 = N*vcov/sigma_squared
	#print(head(vcov))
	hh_inv3 = N*vcov/sigma_squared

	#print("doner")

	#print(hh_inv)
	#sigma_squared = (summary(fit_full)$sigma)^2
	#hh_inv2 = ((N*summary(fit_full)$vcov)/(sigma_squared))
	#print(head(hh_inv))
	#print(head(hh_inv2))
	#print(sqrt(diag(summary(fit_full)$vcov)/sigma_))
	#print(((summary(fit_full)$sigma)))
	#print(sqrt(diag(hh_inv*sigma_squared)))

	##############standard_std_err = sqrt(diag(hh_inv)*0.156959)
	huber_se = get_huber_se(H, hh_inv3, fit_full, expr, N, 1000)


	##############standard_std_err = sqrt(diag(hh_inv)*0.156959)
	#huber_se = get_huber_se(H, hh_inv, fit_full, expr, N, 1000)
	#huber_se_big = get_huber_se(H, hh_inv, fit_full, expr, N, 20000)
	#print(huber_se)

	#lrt <- anova(fit_null,fit_full)

	#aggregate_pvalue <- lrt[[8]][2]
	#print(summary(fit_full))

	# extract coefficients
	coefs <- data.frame(coef(summary(fit_full)))

	#print(coefs)
	#print(tail(coefs))
	coefs$Std..Error = huber_se
	coefs$t.value = coefs$Estimate/coefs$Std..Error
	#print(tail(coefs))
	# use normal distribution to approximate p-value
	coefficient_pvalues <- 2 * (1 - pnorm(abs(coefs$t.value)))

	num_cov = dim(covariates)[2]

	lf_interaction_coefficient_pvalues = coefficient_pvalues[(3+num_cov + num_lf):length(coefficient_pvalues)]
	betas <- coefs[(3+num_cov + num_lf):length(coefficient_pvalues),1]
	std_err <- coefs[(3+num_cov + num_lf):length(coefficient_pvalues),2]

	marginal_beta = coefs[2, 1]
	marginal_std_err = coefs[2,2]
	marginal_pvalue = coefficient_pvalues[2]


	return(list(pvalue=lf_interaction_coefficient_pvalues, beta=betas, std_err=std_err, marginal_beta=marginal_beta, marginal_std_err=marginal_std_err, marginal_pvalue=marginal_pvalue))
}


pass_maf_filter <- function(geno, thresh) {
	booly = TRUE
	num_samples = length(geno)
	af = sum(geno)/(2.0*num_samples)
	if (af > .5) {
		maf = 1.0 - af
	} else {
		maf = af
	}
	if (maf < thresh) {
		booly = FALSE
	}
	return(booly)
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
sample_permutation_file = args[8]
job_number = as.numeric(args[9])
num_jobs = as.numeric(args[10])
total_lines = as.numeric(args[11]) - 1

print(total_lines)


# Determine number of lines each parrallelized job will complete
lines_per_job = ceiling(total_lines/num_jobs)
print(lines_per_job)
start_num = job_number*lines_per_job
end_num = (job_number + 1)*lines_per_job


covariates <- as.matrix(read.table(covariate_file, header=FALSE))
lfs <- as.matrix(read.table(interaction_factor_file, header=FALSE))
groups <- read.table(sample_overlap_file, header=FALSE)$V1
perm <- read.table(sample_permutation_file, header=FALSE)$V1 +1




num_lfs <- dim(lfs)[2]

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
		norm_geno = (geno - mean(geno))/(sd(geno))
		perm_geno = geno[perm]

		norm_perm_geno = (perm_geno - mean(perm_geno))/(sd(perm_geno))

		# Run eqtl analysis
		line_info <- strsplit(line_test,'\t')[[1]]
		ensamble_id = line_info[1]
		rs_id = line_info[2]

		new_line <- paste0(rs_id, "\t", ensamble_id)
		tryCatch(
			{
				#lm_results = run_lf_interaction_eqtl_lm(expr, norm_geno, norm_perm_geno, covariates, lfs, num_lfs)
				lm_results = run_lf_interaction_eqtl_lmm(expr, norm_geno, norm_perm_geno, covariates, lfs, groups, num_lfs)
				new_line <- paste0(new_line, "\t", lm_results$marginal_beta, "\t", lm_results$marginal_std_err, "\t", lm_results$marginal_pvalue)
				for (lf_num in 1:num_lfs) {
					new_line <- paste0(new_line,"\t",lm_results$beta[lf_num], "\t", lm_results$std_err[lf_num], "\t", lm_results$pvalue[lf_num])
				}
				cat(paste0(new_line, "\n"))
        	},
        	error = function(e) {
        		new_line <- paste0(new_line, "\t", 0.0 ,"\t", 1.0, "\t", 1.0)
				for (lf_num in 1:num_lfs) {
					new_line <- paste0(new_line, "\t", 0.0 ,"\t", 1.0, "\t", 1.0)
				}
				cat(paste0(new_line, "\n"))
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