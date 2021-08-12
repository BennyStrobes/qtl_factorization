library(lme4)


y <- read.table("y.txt")
g <- read.table("g.txt")
z <- read.table("z.txt")
cov <- read.table("cov.txt")
y <- y[,1]
g <- g[,1]
cov <- as.matrix(cov)
z <- z[,1]

fit <- lmer(y ~ 0 + g + cov + (1 | z))

convergence_code <- fit@optinfo$conv$opt

if (convergence_code != 0.0) {
	print('LMM DID NOT CONVERGE!')
}


residz = as.numeric(resid(fit))
random_effects_var = as.data.frame(VarCorr(fit))[1,4]
resid_var = as.data.frame(VarCorr(fit))[2,4]

F_mu = coef(summary(fit))[1,1]
C_mu = coef(summary(fit))[2:(dim(cov)[2]+1),1]

alpha_mu = coef(fit)[[1]][,1]


#############
# error checking
# pred = alpha_mu[factor(z)] + g*F_mu + (as.matrix(cov)%*%as.matrix(C_mu))[,1]
# recomputed_resid = y - pred
###############


write.table(resid_var, "residual_variance.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(random_effects_var, "re_variance.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(F_mu, "F_mu.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(C_mu, "C_mu.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(alpha_mu, "alpha_mu.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)