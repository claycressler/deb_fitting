library(car)
## load parameter fits
pars <- vector()
for (geno in c('113-A','2C-AA_take2','322-A','324-A','413-A','4B-AA_take2','711-A','715-A','725-A','815-A','822-A','GLEN-AA_take2','RAMSEY-AA_take2')) {
    name <- load(paste0('Clay_',geno,'_best_fit_parameters.rda'))
    if (name!='p') stop(paste0('parameters called ',name))
    print(sum(!is.na(p[,1])))
    pars <- rbind(pars, cbind(p[!is.na(p[,1]),1:14], genotype=rep(geno,sum(!is.na(p[,1])))))
}

pairs(log10(pars[,1:7]), upper.panel=panel.smooth, lower.panel=NULL, labels=c(expression(alpha),expression(kappa),expression(beta),expression(a),expression(b),expression(c),expression(d)))
dev.copy2pdf(file='Pairwise_parameter_correlations.pdf')

## Life history parameters as a matrix; including initial conditions
## causes the calculation to crash
boot.manova <- function(data, indices) {
    pars <- data[indices,]
    lh.pars <- cbind(pars$alpha,
                     pars$kappa,
                     pars$beta,
                     pars$W.scalar,
                     pars$L.scalar,
                     pars$Rm.scalar,
                     pars$time.scalar)

    ## fit the linear model
    fit.lm <- lm(lh.pars~pars$genotype)
    ## manova
    fit.manova <- Manova(fit.lm)
    ## G matrix
    G <- fit.manova[[1]][[1]]/fit.manova$df
    ## R matrix
    R <- fit.manova$SSPE/fit.manova$error.df
    c(as.vector(cov2cor(G)),as.vector(cov2cor(R)))
}
library(boot)
x <- boot(pars, boot.manova, 10000)

ci.0 <- vector(length=98)
for (i in 1:98) {
    y <- boot.ci(x, type='perc', index=i)
    if (is.null(y)) ci.0[i] <- NA
    else ci.0[i] <- (y[[4]][1,4]/y[[4]][1,5] < 0)
}
matrix(ci.0[1:49], nrow=7, ncol=7)
matrix(ci.0[50:98], nrow=7, ncol=7)
save(x, file='Bootstrapped_G_matrix.R')
