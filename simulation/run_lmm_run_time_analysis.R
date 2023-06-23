args = commandArgs(trailingOnly=TRUE)
library(lme4)
library(lmtest)
library(ggplot2)
library(reshape2)
options(bitmapType='cairo')



create_sample_repeat_structure <- function(sample_size, num_samples_per_individual) {
	lengther = sample_size/num_samples_per_individual
	groups <- rep(1:lengther, each=num_samples_per_individual)
	return(groups)

}

run_lf_interaction_eqtl_lmm <- function(expr, geno, covariates, lfs, groups, num_lf) {
	fit_full <- lmer(expr ~ geno + covariates + lfs + lfs:geno + (1 | groups), REML=FALSE)
	fit_null <- lmer(expr ~ geno + covariates + lfs + (1 | groups), REML=FALSE)

	lrt <- anova(fit_null,fit_full)

	lrt_pvalue <- lrt[[8]][2]
	lrt_fstat <- lrt[[6]][2]
	lrt_df <- lrt[[7]][2]


	# extract coefficients
	coefs <- data.frame(coef(summary(fit_full)))
	# use normal distribution to approximate p-value
	coefficient_pvalues <- 2 * (1 - pnorm(abs(coefs$t.value)))

	num_cov = dim(covariates)[2]

	lf_interaction_coefficient_pvalues = coefficient_pvalues[(3+num_cov + num_lf):length(coefficient_pvalues)]
	betas <- coefs[(3+num_cov + num_lf):length(coefficient_pvalues),1]
	std_err <- coefs[(3+num_cov + num_lf):length(coefficient_pvalues),2]

	marginal_beta = coefs[2, 1]
	marginal_std_err = coefs[2,2]
	marginal_pvalue = coefficient_pvalues[2]


	return(list(pvalue=lf_interaction_coefficient_pvalues, beta=betas, std_err=std_err, marginal_beta=marginal_beta, marginal_std_err=marginal_std_err, marginal_pvalue=marginal_pvalue, lrt_fstat=lrt_fstat, lrt_df=lrt_df, lrt_pvalue=lrt_pvalue))
}


sim_data_dir <- args[1]
eqtl_results_dir <-args[2]
run_time_iter <- args[3]



num_tests=50
simulated_factor = 8
t_statistic=.5
missingness_fraction=.3
num_samples_per_individual=20
n_cov=20
re_var = .05

output_file <- paste0(eqtl_results_dir, "eqtl_lmm_run_time_", run_time_iter, ".txt")
sink(output_file)
new_line <- paste0("sample_size\ttest_iter\trun_time\n")
cat(new_line)

num_samples_input_vector = c(500, 1000, 5000, 10000, 20000)

for (sample_size_iter in 1:length(num_samples_input_vector)) {
	sample_size <- num_samples_input_vector[sample_size_iter]
	groups_sim = create_sample_repeat_structure(sample_size, num_samples_per_individual)


	for (test_iter in 1:num_tests) {
		# Simulate data
		U_sim = matrix(rnorm(n=sample_size*simulated_factor,mean=0, sd=1), ncol=simulated_factor)
		V_sim = rnorm(n=simulated_factor, mean=0, sd=t_statistic)*1.0*(runif(simulated_factor) > (1.0-missingness_fraction))
		cov = matrix(rnorm(n=sample_size*n_cov,mean=0, sd=1), ncol=n_cov)
		cov_effects = rnorm(n=n_cov, mean=0, sd=t_statistic/3)
		af = runif(1,.1,.9)
		geno_a1 = 1.0*(runif(sample_size) > (1.0-af))
		geno_a2 = 1.0*(runif(sample_size) > (1.0-af))
		geno = geno_a1 + geno_a2
		beta = rnorm(n=1,mean=0,sd=.1)
		expr_mean = beta*geno 
		for (factor_iter in 1:simulated_factor) {
			expr_mean = expr_mean + U_sim[,factor_iter]*V_sim[factor_iter]*geno
		}
		for (cov_iter in 1:n_cov) {
			expr_mean = expr_mean + cov[,cov_iter]*cov_effects[cov_iter]
		}
		unique_groups <- unique(groups_sim)
		n_groups = length(unique_groups)
		re_intercepts_sim = rnorm(n=n_groups, mean=0, sd=sqrt(re_var))
		for (group_iter in 1:n_groups) {
			group_name = unique_groups[group_iter]
			group_indices = groups_sim==group_name
			expr_mean[group_indices] = expr_mean[group_indices] + re_intercepts_sim[group_iter]
		}
		expr = rnorm(n=sample_size,mean=expr_mean,sd=1)

		a = Sys.time()
		lm_results = run_lf_interaction_eqtl_lmm(expr, geno, cov, U_sim, groups_sim, dim(U_sim)[2])
		b = Sys.time()
		elapsed_t = paste0(round(as.numeric(difftime(time1 = b, time2 = a, units = "secs")), 3), "")
		new_line=paste0(sample_size, "\t", test_iter, "\t", elapsed_t, "\n")
		cat(new_line)
	}
}
sink()