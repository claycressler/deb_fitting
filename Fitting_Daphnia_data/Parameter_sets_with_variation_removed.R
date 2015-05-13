## Remove the variation in single parameter by fixing at the mean,
## then refit all other parameters. Then calculate the correlation
## between growth and reproduction for the new parameters.

source('trajectory.matching.take2.R')
require(pomp)
require(plyr)

load('Clay_4B-AA_take2_best_fit_parameters.rda')
load('4B-AA_all_growth_egg_data.rda')

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
    np <- array(NA, dim=c(100,15))
    for (i in 1:nrow(p)) {
        if (!is.na(p[i,1])) {
            print(i)
            this.data <- data[[i]]
            this.parset <- p[i,1:14]
            ## Fix the variation
            this.parset[fixed.par] <- median(p[,fixed.par],na.rm=T)
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
                np[i,] <- rep(NA,15)
            else
                np[i,] <- c(x$par, x$lik)
        }
    }
    colnames(np) <- colnames(p)[1:15]

    fixed.parset[[q]] <- np
}
save(fixed.parset, file='Clay_4B-AA_best_fits_w_parameters_fixed.rda')


## Now, I want to calculate the change in the correlation between
## growth and reproduction with these new parameter sets that have had
## their variation in single parameters removed.
## mean von B k value from the original data
mean.k <- mean(p$predK, na.rm=T)

## baseline regressions for partial correlation between K1 and Eggrate
pred.reg.K1.AOD <- lm(p[-20,]$predK1 ~ p[-20,]$predAOD)
pred.reg.Eggrate.AOD <- lm(p[-20,]$predEggrate ~ p[-20,]$predAOD)

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

    reg.K1.AOD <- lm(mod.K1 ~ p$predAOD)
    reg.Eggrate.AOD <- lm(mod.Eggrate ~ p$predAOD)
    new.pcorr.K1.Eggrate <- unname(cor.test(residuals(reg.K1.AOD),residuals(reg.Eggrate.AOD))[[4]])

    ## % change in correlation
    if (base.pcorr.K1.Eggrate < 0)
        delta.pcorr <- -1*(new.pcorr.K1.Eggrate-base.pcorr.K1.Eggrate)/base.pcorr.K1.Eggrate
    else delta.pcorr <- (new.pcorr.K1.Eggrate-base.pcorr.K1.Eggrate)/base.pcorr.K1.Eggrate

    assign(paste0('dCorr.d',par,'.4B'), delta.pcorr)
}

plot.new()
plot.window(xlim=c(0,1), ylim=c(-0.7, 0.7))
at.seq <- seq(0.05,0.95,length=7)
axis(1, at=at.seq, labels=c('alpha','kappa','beta','W sc.', 'L sc.', 'R sc.', 't sc.'))
axis(2)
mtext(side=2, line=3, 'Proportional change in correlation')
box('plot')
abline(h=0)
points(at.seq[1], dCorr.dalpha.4B, pch=21, bg='green')
points(at.seq[2], dCorr.dkappa.4B, pch=21, bg='green')
points(at.seq[3], dCorr.dbeta.4B, pch=21, bg='green')
points(at.seq[4], dCorr.dW.scalar.4B, pch=21, bg='green')
points(at.seq[5], dCorr.dL.scalar.4B, pch=21, bg='green')
points(at.seq[6], dCorr.dRm.scalar.4B, pch=21, bg='green')
points(at.seq[7], dCorr.dtime.scalar.4B, pch=21, bg='green')

dev.copy2pdf(file='Change_in_correlation_from_fixing_parameter.pdf')

