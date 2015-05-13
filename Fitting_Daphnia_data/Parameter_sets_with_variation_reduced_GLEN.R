## Reduce the variation in single parameter by multiplying
## then refit all other parameters. Then calculate the correlation
## between growth and reproduction for the new parameters.

source('trajectory.matching.take2.R')
require(pomp)
require(plyr)

load('Clay_GLEN-AA_take2_best_fit_parameters.rda')
load('GLEN-AA_all_growth_egg_data.rda')

pars.to.fix <- c('alpha','kappa','beta','W.scalar','L.scalar','Rm.scalar','time.scalar')
fixed.parset <- vector(mode='list', length=7)
for (q in 1:7) {
    fixed.par <- pars.to.fix[q]

    fix.pars <- c(c('f','T_M','Rm.0'),fixed.par)
    est.pars <- c('alpha','kappa','beta','W.scalar','L.scalar','Rm.scalar','time.scalar','W.0','L.0','L.sd','R.sd')
    est.pars <- est.pars[-which(est.pars==fixed.par)]
    transform <- c('log','logit',rep('log',12)) ## rules for transforming to the estimation scale

    lower <- c(rep(-50,7),-50,-2,-50,-50)
    upper <- c(rep(50,7),0,0.2,50,50) ## Watch out for the upper bound on W!

    ## Fix a single parameter in turn and refit all others
    np <- array(NA, dim=c(nrow(p),16))
    for (i in 1:nrow(p)) {
        if (!is.na(p[i,1])) {
            print(i)
            this.data <- data[[i]]
            this.parset <- p[i,1:14]
            ## Reduce the variation in the fixed.par by changing this
            ## value by 10% (on the estimation scale, so to change
            ## log(a) by 10%, change a to a^0.9). For kappa, this is a
            ## bit more complicated.
#            if (fixed.par=='kappa') {
#                tK <- 0.95*log(p[i,fixed.par]/(1-p[i,fixed.par]))
#                this.parset[fixed.par] <- exp(tK)/(1+exp(tK))
#            } else this.parset[fixed.par] <- p[i,fixed.par]^0.95
            this.parset[fixed.par] <- p[i,fixed.par]*0.9
            x <- try(trajectory.matching(data=this.data,
                                         pars=this.parset,
                                         transform=transform,
                                         est=est.pars,
                                         method='L-BFGS-B',
                                         lower=lower,
                                         upper=upper,
                                         odefunc="derivs",
                                         initfunc="initmod",
                                         dllname="nondim_deb_take2",
                                         eval.only=FALSE,
                                         maxit=500))
            if (inherits(x, 'try-error'))
                np[i,] <- rep(NA,16)
            else
                np[i,] <- c(x$par, x$lik, x$conv)
        }
    }
    colnames(np) <- c(colnames(p)[1:15],'conv')

    fixed.parset[[q]] <- np
}
save(fixed.parset, file='Clay_GLEN-AA_best_fits_w_variation_reduced_2.rda')


## Now, I want to calculate the change in the correlation between
## growth and reproduction with these new parameter sets that have had
## their variation in single parameters removed.
## mean von B k value from the original data
mean.k <- mean(p$predK, na.rm=T)

## baseline regressions for partial correlation between K1 and Eggrate
pred.reg.K1.AOD <- lm(p[-c(2,50),]$predK1 ~ p[-c(2,50),]$predAOD)
pred.reg.Eggrate.AOD <- lm(p[-c(2,50),]$predEggrate ~ p[-c(2,50),]$predAOD)

## partial correlation coefficients
base.pcorr.K1.Eggrate <- unname(cor.test(residuals(pred.reg.K1.AOD),residuals(pred.reg.Eggrate.AOD))[[4]])

require(deSolve)
dyn.load('nondim_deb_take2.so')

for (j in 1:7) {
    ## Which parameter was fixed?
    par <- pars.to.fix[j]
    ## Store the new Kappa1 and Eggrate values
    mod.K1 <- rep(NA, nrow(p))
    mod.Eggrate <- rep(NA, nrow(p))
    ## Calculate Kappa1 and Eggrate for each of the new parameter sets
    for (i in 1:nrow(p)) {
        if (!is.na(fixed.parset[[j]][i,1])) {
            ## also exclude any parameter sets where conv!=0 {
            if (!(fixed.parset[[j]][i,'conv']!=0)) {
                params <- fixed.parset[[j]][i,1:14]
                y0 <- params[c('W.0','L.0','Rm.0','L.0')]
                names(y0) <- c('W','Lmeas','Rm','L')
                y0[c('Lmeas','L')] <- y0[c('Lmeas','L')]/params['L.scalar']
                y0['W'] <- y0['W']/params['W.scalar']
                times <- c(0,seq(1,p$predAOD[i],2))*params['time.scalar']
                out <- ode(y0,times,func='derivs',dllname='nondim_deb_take2',parms=params,initfunc='initmod')
                out <- as.data.frame(out)
                time <- out$time/params['time.scalar']
                L <- out$Lmeas*params['L.scalar']
                R <- out$Rm-params['beta']
                R[which(R <= 0)] <- 0
                R <- R*params['Rm.scalar']
                if (max(R) > 0) {
                    if (time[min(which(R > 0))]!=max(time))
                        mod.Eggrate[i] <- max(R)/(max(time)-time[min(which(R > 0))]) # ind'l matured on day of death
                    else mod.Eggrate[i] <- max(R)/1.5
                }
                else mod.Eggrate[i] <- 0
                fit <- nls(L ~ Linf*(1-exp(-mean.k*time)), start=c(Linf=3))
                mod.K1[i] <- mean.k*unname(coef(fit)[1])
            }
        }
    }

    ## Exclude datasets that fit the data poorly
    drop <- which(fixed.parset[[j]][,'loglik']-p$loglik > 10)
    reg.K1.AOD <- lm(mod.K1[-drop] ~ p$predAOD[-drop])
    reg.Eggrate.AOD <- lm(mod.Eggrate[-drop] ~ p$predAOD[-drop])
    new.pcorr.K1.Eggrate <- unname(cor.test(residuals(reg.K1.AOD),residuals(reg.Eggrate.AOD))[[4]])

    ## % change in correlation
    if (base.pcorr.K1.Eggrate < 0)
        delta.pcorr <- -1*(new.pcorr.K1.Eggrate-base.pcorr.K1.Eggrate)/base.pcorr.K1.Eggrate
    else delta.pcorr <- (new.pcorr.K1.Eggrate-base.pcorr.K1.Eggrate)/base.pcorr.K1.Eggrate

    assign(paste0('dCorr.d',par,'.GLEN'), delta.pcorr)
}

dCorrdPar <- c(dCorr.dalpha.GLEN, dCorr.dkappa.GLEN, dCorr.dbeta.GLEN, dCorr.dW.scalar.GLEN, dCorr.dL.scalar.GLEN, dCorr.dRm.scalar.GLEN, dCorr.dtime.scalar.GLEN)
names(dCorrdPar) <- c('alpha','kappa','beta','W.scalar','L.scalar','Rm.scalar','time.scalar')

save(dCorrdPar, file='Change_in_correlation_with_reduced_variation_GLEN_2.rda')

