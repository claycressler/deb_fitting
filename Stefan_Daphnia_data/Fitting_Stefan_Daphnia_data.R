for (geno in c('113','322','324','413','711')) {
    ## This new fitting algorithm uses constraints and fits *dimensional*
    ## initial conditions, rather than nondimensional. That is because I
    ## can put reasonable constraints on L.0 and W.0 on the dimensional,
    ## but not nondimensional, scale. I will have to modify the algorithm
    ## to nondimensionalize these initial conditions inside the objective
    ## function. I am also going to remove Lmeas.0 as an initial condition
    ## and simply force it to be equal to L.0.

    ## load all of the Daphnia data for genotype 4B at the highest food level
    load(paste0(geno,'-A_all_growth_egg_data.rda'))

## use trajectory matching to find the best-fit parameter set for each
## individual, beginning from a global search
require(pomp)
require(plyr)

fix.pars <- c('f','T_M','Rm.0')
est.pars <- c('alpha','kappa','beta','W.scalar','L.scalar','Rm.scalar','time.scalar','W.0','L.0','L.sd','R.sd')
transform <- c('log','logit',rep('log',12)) ## rules for transforming to the estimation scale

if (!file.exists(paste0('results-',geno,'-A-stage-one.rda'))) {
    ## create a single set of parameter guesses to try on each individual
    ## to begin the global search.
    box <- cbind(
                 lower=c(alpha=0.5,kappa=0.1,beta=0.01,
                 W.scalar=5e-5,L.scalar=0.1,Rm.scalar=0.1,time.scalar=0.01,
                 f=0.05,T_M=2,
                 W.0=1e-9,Rm.0=1e-9,L.0=0.2,
                 L.sd=0.0001,R.sd=0.0001),
                 upper=c(alpha=50,kappa=0.9,beta=5,
                 W.scalar=5e-3,L.scalar=20,Rm.scalar=20,time.scalar=3,
                 f=0.05,T_M=2,
                 W.0=5,Rm.0=1e-9,L.0=1.2,
                 L.sd=0.1,R.sd=1)
                 )
    set.seed(12340987)
    guesses <- sobolDesign(lower=box[,'lower'],upper=box[,'upper'],nseq=1000)
    guesses$id <- seq_len(nrow(guesses))
    joblist <- dlply(guesses, ~id, unlist)
    ## strip out id
    for (i in 1:length(joblist))
        joblist[[i]] <- joblist[[i]][-which(names(joblist[[i]])=='id')]

    this.data <- data[[1]]
    this.parset <- joblist[[1]]

    source('trajectory.matching.take2.R')
    trajectory.matching(data=this.data, pars=this.parset, transform=transform,
                        est=est.pars, method='subplex',
                        odefunc="derivs", initfunc="initmod",
                        dllname="nondim_deb_take2",
                        eval.only=TRUE)

    ## double check by running outside of trajectory.matching
    y0 <- this.parset[c('W.0','L.0','Rm.0','L.0')]
    names(y0) <- c('W','Lmeas','Rm','L')
    y0['W'] <- y0['W']/this.parset['W.scalar']
    y0[c('Lmeas','L')] <- y0[c('Lmeas','L')]/this.parset['L.scalar']
    times <- this.data$age*this.parset['time.scalar']
    x <- try(ode(y0,c(0,times),func='derivs',parms=this.parset,dllname='nondim_deb_take2',initfunc='initmod'))
    ## remove the first datapoint
    x <- x[-1,]
    ## redimensionalize
    x[,'time'] <- x[,'time']/this.parset['time.scalar']
    x[,3] <- x[,3]*this.parset['L.scalar']
    x[,4] <- x[,4]-this.parset['beta']
    x[which(x[,4]<=0),4] <- 0
    x[,4] <- x[,4]*this.parset['Rm.scalar']
    L.sd <-this.parset['L.sd']
    R.sd <-this.parset['R.sd']
    lik <- -sum(dnorm(this.data$Lobs, x[,3], L.sd, log=T) + dnorm(this.data$Robs, x[,4], R.sd, log=T), na.rm=T)

    results = vector(mode='list', length=length(data))
    for (i in 1:length(data)) {
        this.data <- data[[i]]
        if (nrow(this.data) > 7) {
            lik <- lapply(joblist,
                          function(x) trajectory.matching(data=this.data,
                                                          pars=x,
                                                          transform=transform,
                                                          est=est.pars,
                                                          odefunc="derivs",
                                                          initfunc="initmod",
                                                          dllname="nondim_deb_take2",
                                                          eval.only=TRUE))
            this.results <- ldply(lik, function(x) c(x$par, loglik=x$lik))
            this.results <- this.results[order(this.results$loglik),]
            results[[i]] <- this.results
        }
    }
    save(results, file=paste0('results-',geno,'-A-stage-one.rda'))
} else load(paste0('results-',geno,'-A-stage-one.rda'))

## Since I am going to be using optimization with constraints, I need
## to specify the upper and lower bounds for searching.
lower <- c(rep(-50,7),-50,-2,-50,-50)
upper <- c(rep(50,7),2,0.2,50,50)

## Take the top 100 parameter sets for each individual and do a full
## trajectory.matching run, using L-BFGS-B with only 500 steps initially.
results2 <- vector(mode='list', length=length(results))
for (i in 1:length(results)) {
    if (is.data.frame(results[[i]])) {
        this.data <- data[[i]]
        vals <- array(NA, dim=c(100, 17))
        for (j in 1:100) {
            print(paste('dataset',i,'parameter set',j))
            this.parset <- results[[i]][j,2:15] ## extract the parameter values
            start.loglik <- results[[i]][j,16]
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
            if (inherits(x,'try-error'))
                vals[j,] <- c(as.numeric(this.parset), loglik=NA, conv=NA, dloglik=NA)
            else {
                finish.loglik <- x$lik
                loglik.diff <- start.loglik-finish.loglik
                vals[j,] <- c(x$par, loglik=x$lik, conv=x$conv, dloglik=loglik.diff)
            }
        }
        vals <- as.data.frame(vals)
        colnames(vals) <- c('alpha','kappa','beta','W.scalar','L.scalar','Rm.scalar','time.scalar','f','T_M','W.0','Rm.0','L.0','L.sd','R.sd','loglik','conv','dloglik')
        results2[[i]] <- vals
        save(results2, file=paste0('results-',geno,'-A-stage-two.rda'))
    }
}

## looking at results2 - individuals 47 and 67 couldn't be fit; both of these are long-lived individuals that never reproduced - skip these. Skip all individuals whose minimum likelihood is either Inf or -Inf
inds <- which(is.finite(sapply(results2, function(x) min(x$loglik, na.rm=T))))
results3 <- vector(mode='list', length=length(results2))
for (i in inds) {
  this.data <- data[[i]]
  vals <- array(NA, dim=c(100, 17))
  for (j in 1:100) {
      print(paste('dataset', i, 'parset', j))
      this.parset <- results2[[i]][j,1:14] ## extract the parameter values
      start.loglik <- as.numeric(results2[[i]][j,15])

      if (is.na(start.loglik))
          vals[j,] <- c(as.numeric(this.parset), loglik=start.loglik, conv=NA, dloglik=NA)

      else {
          count <- 1
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

          if (inherits(x,'try-error'))
              vals[j,] <- c(as.numeric(this.parset), loglik=start.loglik, conv=NA, dloglik=0)
          else if (is.na(x$lik))
              vals[j,] <- c(as.numeric(this.parset), loglik=start.loglik, conv=NA, dloglik=0)
          else {
              finish.loglik <- x$lik
              loglik.diff <- start.loglik-finish.loglik
              vals[j,] <- c(as.numeric(x$par), loglik=x$lik, conv=x$conv, dloglik=loglik.diff)
              while (abs(loglik.diff) > 0.1) {
                  count <- count+1
                  print(paste('dloglik', loglik.diff, 'extra TM steps', 500*(count-1)))
                  start.loglik <- finish.loglik
                  this.parset <- x$par
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
                  if (inherits(x,'try-error')) {
                      vals[j,] <- c(as.numeric(this.parset), loglik=start.loglik, conv=NA, dloglik=loglik.diff)
                      loglik.diff <- 0
                  }
                  else if (is.na(x$lik)) {
                      vals[j,] <- c(as.numeric(this.parset), loglik=start.loglik, conv=NA, dloglik=loglik.diff)
                      loglik.diff <- 0
                  }
                  else {
                      finish.loglik <- x$lik
                      loglik.diff <- start.loglik-finish.loglik
                      vals[j,] <- c(as.numeric(x$par), loglik=x$lik, conv=x$conv, dloglik=loglik.diff)
                  }
              }
          }
      }
  }
  vals <- as.data.frame(vals)
  colnames(vals) <- c('alpha','kappa','beta','W.scalar','L.scalar','Rm.scalar','time.scalar','f','T_M','W.0','Rm.0','L.0','L.sd','R.sd','loglik','conv','dloglik')
  results3[[i]] <- vals
  save(results3, file=paste0('results-',geno,'-A-stage-three.rda'))
}
}



