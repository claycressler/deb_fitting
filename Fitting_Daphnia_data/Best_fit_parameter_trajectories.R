## Plot trajectories for each parameter set as well as the real
## data. That allows me to compare the fits to the truth, and throw
## out any parameter sets that fit shitty.
require(deSolve)
dyn.load('nondim_deb_take2.so')

for (geno in c('2C','4B','GLEN','RAMSEY')) {
    load(paste0('Clay_',geno,'-AA_take2_best_fit_parameters.rda'))
    load(paste0(geno,'-AA_all_growth_egg_data.rda'))
    for (i in 1:nrow(p)) {
        if (!is.na(p[i,1])) {# check to make sure the data could be fit
            y0 <- p[i,c('W.0','L.0','Rm.0','L.0')]
            names(y0) <- c('W','Lmeas','Rm','L')
            y0[c('Lmeas','L')] <- y0[c('Lmeas','L')]/p[i,'L.scalar']
            y0['W'] <- y0['W']/p[i,'W.scalar']
            times <- c(0,data[[i]]$age*p[i,'time.scalar'])
            params <- p[i,1:14]
            out <- ode(y0,times,func='derivs',dllname='nondim_deb_take2',parms=params,initfunc='initmod')
            out <- as.data.frame(out)
            time <- out$time/p[i,'time.scalar']
            L <- out$Lmeas*p[i,'L.scalar']
            R <- out$Rm-p[i,'beta']
            R[which(R <= 0)] <- 0
            R <- R*p[i,'Rm.scalar']

            par(mfrow=c(1,2), oma=rep(0.5,4), mar=c(5,5,2,0))
            plot(data[[i]]$age, data[[i]]$Lobs, xlim=c(0,max(data[[i]]$age)), pch=21, bg=1, main=paste0(geno,' indl ',i), xlab='Age',ylab='Length')
            lines(time, L, col=2, lwd=1.5)
            plot(data[[i]]$age, data[[i]]$Robs, xlim=c(0,max(data[[i]]$age)), pch=21, bg=1, ylim=c(0,max(c(data[[i]]$Robs,R),na.rm=T)),xlab='Age',ylab='Cum. Eggs')
            lines(time, R, col=2, lwd=1.5)
            dev.copy2pdf(file=paste0('Data_and_fits/',geno,'-AA_indl_',i,'.pdf'))
        }
    }
}

## Looking at the plots, I can safely remove a number of poor-fitting animals.
