## Load the best fitting parameter for each 2C-AA individual and load
## the observed data.
load('Clay_2C-AA_best_fit_parameters.rda')
load('2C-AA_all_growth_egg_data.rda')

require(deSolve)
dyn.load('nondim_deb.so')

for (i in 1:nrow(p.2C)) {
    if (!is.na(p.2C[i,1])) {
        y0 <- p.2C[i,c('W.0','Lmeas.0','Rm.0','L.0')]
        times <- c(0,data[[i]]$age*p.2C[i,'time.scalar'])
        params <- p.2C[i,1:15]
        out <- ode(y0,times,func='derivs',dllname='nondim_deb',parms=params,initfunc='initmod')
        out <- as.data.frame(out)
        ## re-dimensionalize
        time <- out$time/p.2C[i,'time.scalar']
        L <- out$Lmeas.0*p.2C[i,'L.scalar']
        R <- out$Rm.0-p.2C[i,'beta']
        R[which(R <= 0)] <- 0
        R <- R*p.2C[i,'Rm.scalar']
        params['W.0'] <- params['W.0']*params['W.scalar']
        params['L.0'] <- params['L.0']*params['L.scalar']
        params['Lmeas.0'] <- params['Lmeas.0']*params['L.scalar']

        par(mfrow=c(1,2), oma=rep(0.5,4), mar=c(5,5,2,0))
        plot(data[[i]]$age, data[[i]]$Lobs, xlim=c(0,max(data[[i]]$age)), pch=21, bg=1, main=paste('Individual',i), xlab='Age',ylab='Length')
        lines(time, L, col=2, lwd=1.5)
        legend(x='bottomright', paste0(names(params),'=',signif(params,3)),bty='n')

        plot(data[[i]]$age, data[[i]]$Robs, xlim=c(0,max(data[[i]]$age)), pch=21, bg=1, ylim=c(0,max(c(data[[i]]$Robs,R))),xlab='Age',ylab='Cum. Eggs')
        lines(time, R, col=2, lwd=1.5)
        dev.copy2pdf(file=paste0('2C-AA_data+fits/Individual-',i,'.pdf'))
    }
}

