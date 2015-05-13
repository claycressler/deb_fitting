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

    ## Drop parameter sets that fit the data very poorly after fixing
    drop <- which(fixed.parset[[j]][1:nrow(p),'loglik']-p$loglik > 20)

    reg.K1.AOD <- lm(mod.K1[-drop] ~ p$predAOD[-drop])
    reg.Eggrate.AOD <- lm(mod.Eggrate[-drop] ~ p$predAOD[-drop])
    new.pcorr.K1.Eggrate <- unname(cor.test(residuals(reg.K1.AOD),residuals(reg.Eggrate.AOD))[[4]])

    ## % change in correlation
    if (base.pcorr.K1.Eggrate < 0)
        delta.pcorr <- -1*(new.pcorr.K1.Eggrate-base.pcorr.K1.Eggrate)/base.pcorr.K1.Eggrate
    else delta.pcorr <- (new.pcorr.K1.Eggrate-base.pcorr.K1.Eggrate)/base.pcorr.K1.Eggrate

    assign(paste0('dCorr.d',par,'.4B'), delta.pcorr)
}

dCorrdPar <- c(dCorr.dalpha.4B, dCorr.dkappa.4B, dCorr.dbeta.4B, dCorr.dW.scalar.4B, dCorr.dL.scalar.4B, dCorr.dRm.scalar.4B, dCorr.dtime.scalar.4B)
names(dCorrdPar) <- c('alpha','kappa','beta','W.scalar','L.scalar','Rm.scalar','time.scalar')

save(dCorrdPar, file='Change_in_correlation_with_fixing_parameters_4B.rda')

load('Change_in_correlation_with_fixing_parameters_2C.rda')
dCorrdPar.2C <- dCorrdPar
load('Change_in_correlation_with_fixing_parameters_4B.rda')
dCorrdPar.4B <- dCorrdPar
load('Change_in_correlation_with_fixing_parameters_GLEN.rda')
dCorrdPar.GLEN <- dCorrdPar
load('Change_in_correlation_with_fixing_parameters_RAMSEY.rda')
dCorrdPar.RAMSEY <- dCorrdPar

par(oma=c(4,4,2,2), mar=rep(0,4))
plot.new()
plot.window(xlim=c(0,1), ylim=c(-3, 2))
at.seq <- seq(0.05,0.95,length=7)
axis(1, at=at.seq, labels=c('alpha','kappa','beta','W sc.', 'L sc.', 'R sc.', 't sc.'))
axis(2)
mtext(side=2, line=3, 'Proportional change in correlation')
box('plot')
abline(h=0)
points(at.seq[1], dCorrdPar.4B['alpha'], pch=21, bg='green')
points(at.seq[1], dCorrdPar.2C['alpha'], pch=21, bg='blue')
points(at.seq[1], dCorrdPar.GLEN['alpha'], pch=21, bg='gold')
points(at.seq[1], dCorrdPar.RAMSEY['alpha'], pch=21, bg='red')

points(at.seq[2], dCorrdPar.4B['kappa'], pch=21, bg='green')
points(at.seq[2], dCorrdPar.2C['kappa'], pch=21, bg='blue')
points(at.seq[2], dCorrdPar.GLEN['kappa'], pch=21, bg='gold')
points(at.seq[2], dCorrdPar.RAMSEY['kappa'], pch=21, bg='red')

points(at.seq[3], dCorrdPar.4B['beta'], pch=21, bg='green')
points(at.seq[3], dCorrdPar.2C['beta'], pch=21, bg='blue')
points(at.seq[3], dCorrdPar.GLEN['beta'], pch=21, bg='gold')
points(at.seq[3], dCorrdPar.RAMSEY['beta'], pch=21, bg='red')

points(at.seq[4], dCorrdPar.4B['W.scalar'], pch=21, bg='green')
points(at.seq[4], dCorrdPar.2C['W.scalar'], pch=21, bg='blue')
points(at.seq[4], dCorrdPar.GLEN['W.scalar'], pch=21, bg='gold')
points(at.seq[4], dCorrdPar.RAMSEY['W.scalar'], pch=21, bg='red')

points(at.seq[5], dCorrdPar.4B['L.scalar'], pch=21, bg='green')
points(at.seq[5], dCorrdPar.2C['L.scalar'], pch=21, bg='blue')
points(at.seq[5], dCorrdPar.GLEN['L.scalar'], pch=21, bg='gold')
points(at.seq[5], dCorrdPar.RAMSEY['L.scalar'], pch=21, bg='red')

points(at.seq[6], dCorrdPar.4B['Rm.scalar'], pch=21, bg='green')
points(at.seq[6], dCorrdPar.2C['Rm.scalar'], pch=21, bg='blue')
points(at.seq[6], dCorrdPar.GLEN['Rm.scalar'], pch=21, bg='gold')
points(at.seq[6], dCorrdPar.RAMSEY['Rm.scalar'], pch=21, bg='red')

points(at.seq[7], dCorrdPar.4B['time.scalar'], pch=21, bg='green')
points(at.seq[7], dCorrdPar.2C['time.scalar'], pch=21, bg='blue')
points(at.seq[7], dCorrdPar.GLEN['time.scalar'], pch=21, bg='gold')
points(at.seq[7], dCorrdPar.RAMSEY['time.scalar'], pch=21, bg='red')

legend(x='bottomright', legend=c('4B','2C','GLEN','RAMSEY'), pch=21, pt.bg=c('green','blue','gold','red'), bty='n')

dev.copy2pdf(file='Change_in_correlation_from_fixing_parameter.pdf')

