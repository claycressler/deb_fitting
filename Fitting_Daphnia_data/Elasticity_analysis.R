for (geno in c('2C','4B','GLEN','RAMSEY')) {

    load(paste0('Clay_',geno,'-AA_take2_best_fit_parameters.rda'))
    ## Throw away some bad-fitting parameter sets
    if (geno=='4B') p <- p[-20,]
    else if (geno=='GLEN') p <- p[-c(2,50),]

    ## What is the mean(k) value? Needed for re-fitting Linf to the
    ## modified parameter sets
    mean.k <- mean(p$predK,na.rm=T)

    ## baseline regressions for partial correlation between K1 and Eggrate
    pred.reg.K1.AOD <- lm(p$predK1 ~ p$predAOD)
    pred.reg.Eggrate.AOD <- lm(p$predEggrate ~ p$predAOD)

    ## partial correlation coefficients
    base.pcorr.K1.Eggrate <- unname(cor.test(residuals(pred.reg.K1.AOD),residuals(pred.reg.Eggrate.AOD))[[4]])

    require(deSolve)
    dyn.load('nondim_deb_take2.so')
    ## Vary parameters of interest (alpha, kappa, beta, a, b, c, d) by
    ## a small amount, simulate a new trajectory, recompute Linf and
    ## Eggrate, and then recalculate the correlation
    for (par in c('alpha','kappa','beta','W.scalar','L.scalar','Rm.scalar','time.scalar')) {
        mod.Linf <- rep(NA, nrow(p))
        mod.K1 <- rep(NA, nrow(p))
        mod.Eggrate <- rep(NA, nrow(p))
        for (i in 1:nrow(p)) {
            if (!is.na(p[i,1])) {
                params <- as.numeric(p[i,1:14])
                names(params) <- colnames(p)[1:14]
                params[par] <- params[par]-0.01*params[par]
                y0 <- as.numeric(params[c('W.0','L.0','Rm.0','L.0')])
                names(y0) <- c('W','Lmeas','Rm','L')
                y0[c('Lmeas','L')] <- y0[c('Lmeas','L')]/as.numeric(params['L.scalar'])
                y0['W'] <- y0['W']/as.numeric(params['W.scalar'])
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
                mod.Linf[i] <- unname(coef(fit)[1])
                mod.K1[i] <- mean.k*unname(coef(fit)[1])
            }
        }
        ## Compute the new correlation coefficient between K1 and Eggrate
        reg.K1.AOD <- lm(mod.K1 ~ p$predAOD)
        reg.Eggrate.AOD <- lm(mod.Eggrate ~ p$predAOD)
        new.pcorr.K1.Eggrate <- unname(cor.test(residuals(reg.K1.AOD),residuals(reg.Eggrate.AOD))[[4]])

        ## % change in correlation
        if (base.pcorr.K1.Eggrate < 0)
            delta.pcorr <- -1*(new.pcorr.K1.Eggrate-base.pcorr.K1.Eggrate)/base.pcorr.K1.Eggrate
        else delta.pcorr <- (new.pcorr.K1.Eggrate-base.pcorr.K1.Eggrate)/base.pcorr.K1.Eggrate

        assign(paste0('dCorr.d',par,'.',geno), delta.pcorr)
    }
}

for (geno in c('113','322','324','715','725','815','822')) {
    load(paste0('Clay_',geno,'-A_best_fit_parameters.rda'))

    ## What is the mean(k) value? Needed for re-fitting Linf to the
    ## modified parameter sets
    mean.k <- mean(p$predK,na.rm=T)

    ## baseline regressions for partial correlation between K1 and Eggrate
    pred.reg.K1.AOD <- lm(p$predK1 ~ p$predAOD)
    pred.reg.Eggrate.AOD <- lm(p$predEggrate ~ p$predAOD)

    ## partial correlation coefficients
    base.pcorr.K1.Eggrate <- unname(cor.test(residuals(pred.reg.K1.AOD),residuals(pred.reg.Eggrate.AOD))[[4]])

    require(deSolve)
    dyn.load('nondim_deb_take2.so')
    ## Vary parameters of interest (alpha, kappa, beta, a, b, c, d) by
    ## a small amount, simulate a new trajectory, recompute Linf and
    ## Eggrate, and then recalculate the correlation
    for (par in c('alpha','kappa','beta','W.scalar','L.scalar','Rm.scalar','time.scalar')) {
        mod.Linf <- rep(NA, nrow(p))
        mod.K1 <- rep(NA, nrow(p))
        mod.Eggrate <- rep(NA, nrow(p))
        for (i in 1:nrow(p)) {
            if (!is.na(p[i,1])) {
                params <- as.numeric(p[i,1:14])
                names(params) <- colnames(p)[1:14]
                params[par] <- params[par]-0.01*params[par]
                y0 <- as.numeric(params[c('W.0','L.0','Rm.0','L.0')])
                names(y0) <- c('W','Lmeas','Rm','L')
                y0[c('Lmeas','L')] <- y0[c('Lmeas','L')]/as.numeric(params['L.scalar'])
                y0['W'] <- y0['W']/as.numeric(params['W.scalar'])
                times <- c(0,seq(1,p$predAOD[i],2))*params['time.scalar']
                out <- ode(y0,times,func='derivs',dllname='nondim_deb_take2',parms=params,initfunc='initmod')
                out <- as.data.frame(out)
                time <- out$time/p[i,'time.scalar']
                L <- out$Lmeas*p[i,'L.scalar']
                R <- out$Rm-p[i,'beta']
                R[which(R <= 0)] <- 0
                R <- R*p[i,'Rm.scalar']
                if (max(R) > 0) {
                    if (time[min(which(R > 0))]!=max(time))
                        mod.Eggrate[i] <- max(R)/(max(time)-time[min(which(R > 0))]) # ind'l matured on day of death
                    else mod.Eggrate[i] <- max(R)/1.5
                }
                else mod.Eggrate[i] <- 0
                fit <- nls(L ~ Linf*(1-exp(-mean.k*time)), start=c(Linf=3))
                mod.Linf[i] <- unname(coef(fit)[1])
                mod.K1[i] <- mean.k*unname(coef(fit)[1])
            }
        }
        ## Compute the new correlation coefficient between K1 and Eggrate
        reg.K1.AOD <- lm(mod.K1 ~ p$predAOD)
        reg.Eggrate.AOD <- lm(mod.Eggrate ~ p$predAOD)
        new.pcorr.K1.Eggrate <- unname(cor.test(residuals(reg.K1.AOD),residuals(reg.Eggrate.AOD))[[4]])

        ## % change in correlation
        if (base.pcorr.K1.Eggrate < 0)
            delta.pcorr <- -1*(new.pcorr.K1.Eggrate-base.pcorr.K1.Eggrate)/base.pcorr.K1.Eggrate
        else delta.pcorr <- (new.pcorr.K1.Eggrate-base.pcorr.K1.Eggrate)/base.pcorr.K1.Eggrate

        assign(paste0('dCorr.d',par,'.',geno), delta.pcorr)
    }
}

geno.seq <- c('2C','4B','GLEN','RAMSEY','113','322','324','715','725','815','822')
require(RColorBrewer)
col.seq <- brewer.pal(length(geno.seq), 'Set3')

plot.new()
plot.window(xlim=c(0,1),ylim=c(-1.5,0.2))
at.seq <- seq(0.05,0.95,length=7)
axis(1, at=at.seq, labels=c('alpha','kappa','beta','W sc.', 'L sc.', 'R sc.', 't sc.'))
axis(2)
mtext(side=2, line=3, 'Proportional change in correlation')
box('plot')
abline(h=0)
for (g in 1:length(geno.seq))
    points(at.seq[1], get(paste0('dCorr.dalpha.',geno.seq[g])), pch=21, bg=col.seq[g], cex=1.8)

for (g in 1:length(geno.seq))
    points(at.seq[2], get(paste0('dCorr.dkappa.',geno.seq[g])), pch=21, bg=col.seq[g], cex=1.8)

for (g in 1:length(geno.seq))
    points(at.seq[3], get(paste0('dCorr.dbeta.',geno.seq[g])), pch=21, bg=col.seq[g], cex=1.8)

for (g in 1:length(geno.seq))
    points(at.seq[4], get(paste0('dCorr.dW.scalar.',geno.seq[g])), pch=21, bg=col.seq[g], cex=1.8)

for (g in 1:length(geno.seq))
    points(at.seq[5], get(paste0('dCorr.dL.scalar.',geno.seq[g])), pch=21, bg=col.seq[g], cex=1.8)

for (g in 1:length(geno.seq))
    points(at.seq[6], get(paste0('dCorr.dRm.scalar.',geno.seq[g])), pch=21, bg=col.seq[g], cex=1.8)

for (g in 1:length(geno.seq))
    points(at.seq[7], get(paste0('dCorr.dtime.scalar.',geno.seq[g])), pch=21, bg=col.seq[g], cex=1.8)

legend(x='bottomright', legend=geno.seq, pch=21, pt.bg=col.seq, col=col.seq, bty='n')

dev.copy2pdf(file='Elasticity_analysis.pdf')


