for (geno in c('2C','4B','GLEN','RAMSEY')) {
    load(paste0('results-',geno,'-AA-take2-stage-three.rda'))
    p <- array(NA, dim=c(length(results3),15))
    for (i in 1:length(results3)) {
        if (length(results3[[i]])>0)
            p[i,] <- as.numeric(results3[[i]][order(results3[[i]][,'loglik']),][1,1:15])
    }
    colnames(p) <- c('alpha','kappa','beta','W.scalar','L.scalar','Rm.scalar','time.scalar','f','T_M','W.0','Rm.0','L.0','L.sd','R.sd','loglik')
    save(p, file=paste0('Clay_',geno,'-AA_take2_best_fit_parameters.rda'))
}

## Add the Linf value Adriana computed using NLS, the observed
## eggrate, and the observed age at death to the list of best fit
## parameter values
x <- read.csv('~/Dropbox/Cressler/DEB_fitting/Adriana_Daphnia_data/IndivData - WITH linf.csv')
y <- read.csv('~/Dropbox/Cressler/DEB_fitting/Adriana_Daphnia_data/IndivData - NO linf.csv')

for (geno in c('2C','4B','GLEN','RAMSEY')) {
    size = read.csv(paste0('~/Dropbox/Cressler/DEB_fitting/Adriana_Daphnia_data/',geno,'-AA-GROWTH.csv'), header=F, na.strings='na', colClasses='character')
    load(paste0('Clay_',geno,'-AA_take2_best_fit_parameters.rda'))

    ## match ID in 'x' with ID in 'size' and then with row in 'p'
    IDrow <- which(size[,1]=='ID')
    Linf <- rep(NA, nrow(p))
    aod <- rep(NA, nrow(p))
    eggrate <- rep(NA, nrow(p))
    for (i in 2:ncol(size)) {
        tID <- size[IDrow,i]
        rx <- which(x$ID==tID & x$Geno==geno & x$Treat==0.025)
        ry <- which(y$ID==tID & y$Geno==geno & y$Treat==0.025)
        if (length(rx) > 0 & length(ry) > 0) { ## Linf and Eggrate were computed for this individual
            Linf[i-1] <- x$Linf[rx]
            aod[i-1] <- x$AOD[rx]
            eggrate[i-1] <- y$Eggrate[ry]
        }
    }
    p <- as.data.frame(p)
    p$Linf <- Linf
    p$AOD <- aod
    p$Eggrate <- eggrate

    save(p, file=paste0('Clay_',geno,'-AA_take2_best_fit_parameters.rda'))
}

## Compute Linf, k, maximum length, mean daily egg rate, and age at
## death from trajectories simulated using the best fit parameter
## values (age at death is added because Adriana's data was missing
## AOD for some individuals, and I don't know why).
require(deSolve)
dyn.load('nondim_deb_take2.so')

for (geno in c('2C','4B','GLEN','RAMSEY')) {
    load(paste0('Clay_',geno,'-AA_take2_best_fit_parameters.rda'))
    load(paste0(geno,'-AA_all_growth_egg_data.rda'))
    predLmax <- rep(NA, nrow(p))
    predLinf <- rep(NA, nrow(p))
    predK <- rep(NA, nrow(p))
    predEggrate <- rep(NA, nrow(p))
    predAOD <- rep(NA, nrow(p))
    for (i in 1:nrow(p)) {
        if(!is.na(p[i,1])) {
            y0 <- as.numeric(p[i,c('W.0','L.0','Rm.0','L.0')])
            names(y0) <- c('W','Lmeas','Rm','L')
            y0[c('Lmeas','L')] <- y0[c('Lmeas','L')]/p[i,'L.scalar']
            y0['W'] <- y0['W']/p[i,'W.scalar']
            times <- c(0,data[[i]]$age*p[i,'time.scalar'])
            params <- as.numeric(p[i,1:14])
            out <- ode(y0,times,func='derivs',dllname='nondim_deb_take2',parms=params,initfunc='initmod')
            out <- as.data.frame(out)
            time <- out$time/p[i,'time.scalar']
            L <- out$Lmeas*p[i,'L.scalar']
            R <- out$Rm-p[i,'beta']
            R[which(R <= 0)] <- 0
            R <- R*p[i,'Rm.scalar']
            predLmax[i] <- max(L) ## maximum length
            if (max(R) > 0)
                predEggrate[i] <- max(R)/(max(time)-time[min(which(R > 0))]) ## mean daily egg rate
            else predEggrate[i] <- 0
            ## Fit von Bertalanffy equation and save fitted Linf and k values
            fit <- nls(L ~ Linf*(1-exp(-k*time)), start=c(Linf=3,k=0.1))
            predLinf[i] <- unname(coef(fit)[1])
            predK[i] <- unname(coef(fit)[2])
            predAOD[i] <- max(time)
        }
    }
    ## add these to the data frame containing the best-fit parameter values
    p <- as.data.frame(p)
    p$predLmax <- predLmax
    p$predEggrate <- predEggrate
    p$predLinf <- predLinf
    p$predK <- predK
    p$predAOD <- predAOD
    save(p, file=paste0('Clay_',geno,'-AA_take2_best_fit_parameters.rda'))
}

## Adriana computed an aggregate parameter K1 = Linf*mean(k), where
## Linf is the NLS fit when k was held fixed at the mean k value from
## independent fitting. For her data, mean(k) = 0.161.
for (geno in c('2C','4B','GLEN','RAMSEY')) {
    load(paste0('Clay_',geno,'-AA_take2_best_fit_parameters.rda'))
    load(paste0(geno,'-AA_all_growth_egg_data.rda'))
    p$K1 <- p$Linf*0.161 ## Adriana's K1 calculation
    predK1 <- rep(NA, nrow(p))
    mean.K <- mean(p$predK,na.rm=T)
    for (i in 1:nrow(p)) {
        if(!is.na(p[i,1])) {
            y0 <- as.numeric(p[i,c('W.0','L.0','Rm.0','L.0')])
            names(y0) <- c('W','Lmeas','Rm','L')
            y0[c('Lmeas','L')] <- y0[c('Lmeas','L')]/p[i,'L.scalar']
            y0['W'] <- y0['W']/p[i,'W.scalar']
            times <- c(0,data[[i]]$age*p[i,'time.scalar'])
            params <- as.numeric(p[i,1:14])
            out <- ode(y0,times,func='derivs',dllname='nondim_deb_take2',parms=params,initfunc='initmod')
            out <- as.data.frame(out)
            time <- out$time/p[i,'time.scalar']
            L <- out$Lmeas*p[i,'L.scalar']
            R <- out$Rm-p[i,'beta']
            R[which(R <= 0)] <- 0
            R <- R*p[i,'Rm.scalar']
            ## Fit von Bertalanffy equation
            fit <- nls(L ~ Linf*(1-exp(-mean.K*time)), start=c(Linf=3))
            predK1[i] <- mean.K*unname(coef(fit)[1])
        }
    }
    ## add these to the data frame containing the best-fit parameter values
    p$predK1 <- predK1
    save(p, file=paste0('Clay_',geno,'-AA_take2_best_fit_parameters.rda'))
}
