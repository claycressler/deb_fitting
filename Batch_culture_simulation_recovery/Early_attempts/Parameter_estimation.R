require(pomp)
require(plyr)

## Code for trajectory matching using the batch transfer ODE model
source('trajectory.matching.R')
dyn.load('my_deb_batch_constantK.so')

## Four different true parameter sets
parameter.order <- c('kappa','km','eG','eR','nu','Rmbar','pam','Fh','eA','vol','E.0','L.0','Re.0','R.0','F.0','L.sd')

if(!file.exists('true_parameters.rda')) {
  true.parameters <- vector(mode='list',length=4)
  true.parameters[[1]] <- c(kappa=0.6, km=0.33, eG=0.0017, eR=0.00868, nu=18.1, Rmbar=0.0189, pam=0.005187, Fh=3.09e-5, eA=0.7, vol=20, E.0=0.005187/18.1*0.85^3, L.0=0.85, Re.0=0, R.0=0, F.0=0.025, L.sd=0.02)
  true.parameters[[2]] <- c(kappa=0.7, km=0.23, eG=0.0017, eR=0.00868, nu=3.1, Rmbar=0.0189, pam=0.005187, Fh=3.09e-5, eA=0.6, vol=20, E.0=0.005187/3.1*0.75^3, L.0=0.75, Re.0=0, R.0=0, F.0=0.025, L.sd=0.02)
  true.parameters[[3]] <- c(kappa=0.4, km=0.1, eG=0.001, eR=0.0048, nu=10, Rmbar=0.0189, pam=0.005187, Fh=3.09e-5, eA=0.5, vol=20, E.0=0.005187/10*0.5^3, L.0=0.5, Re.0=0, R.0=0, F.0=0.05, L.sd=0.02)
  true.parameters[[4]] <- c(kappa=0.5, km=0.15, eG=0.001, eR=0.005, nu=5, Rmbar=0.0189, pam=0.005187, Fh=3.09e-5, eA=1, vol=20, E.0=0.005187/5*1^3, L.0=1, Re.0=0, R.0=0, F.0=0.05, L.sd=0.02)
  save(true.parameters, file='true_parameters.rda')
} else load('true_parameters.rda')

cummin <- function(x) sapply(2:length(x), function(y) x[y]-x[y-1])
if(!file.exists('observed_data.rda')) {
  ## Generate four datasets with measurement error
  observed.data <- vector(mode='list', length=4)
  for (i in 1:4) {
    y0 <- true.parameters[[i]][grep('.0',parameter.order)]
    ## simulate with two-day transfers (specified in 'events')
    transfers <- seq(1,75,2)
    out <- ode(y0, times=seq(0,75), func="derivs", parms=true.parameters[[i]], dllname="my_deb_batch_constantK", initfunc="initmod", events=list(func="event",time=transfers))
    ## simulate measurements
    Lobs <- rnorm(length(out[(transfers+1),3]),out[(transfers+1),3],0.02)
    ## Make sure Robs is non-decreasing
    Robs <- rpois(length(out[(transfers+1),5]),out[(transfers+1),5])
    observed.data[[i]] <- as.data.frame(cbind(time=transfers,cbind(Lobs,Robs)))
  }
  save(observed.data, file='observed_data.rda')
} else load('observed_data.rda')

## Performing a global search of parameter space for each dataset,
## estimating the following parameters:
estimated.pars <- c('kappa','km','eG','eR','nu','Rmbar','pam','Fh','eA','E.0','L.0')

for (i in 1:4) {
## Encompass our initial ignorance of the estimated parameter values by a box
    box <- cbind(
                 lower=c(kappa=0.1, km=0.01, eG=0.0001, eR=0.0001, nu=1, Rmbar=0.001, pam=0.0005, Fh=8e-5, eA=0.5, E.0=0.00001, L.0=0.5),
                 upper=c(kappa=0.9, km=1, eG=0.01, eR=0.01, nu=25, Rmbar=0.1, pam=0.05, Fh=8e-4, eA=1, E.0=0.0001, L.0=1.2)
                 )
    ## Create a Sobol low discrepancy sequence of length 10000
    guesses <- sobolDesign(lower=box[,'lower'],upper=box[,'upper'],nseq=10000)
    guesses$id <- seq_len(nrow(guesses))
    joblist <- dlply(guesses, ~id, unlist)
    ## fixed parameters will be constant
    fixed.pars <- c(vol=20, Re.0=0, R.0=0, F.0=true.parameters[[i]]['F.0'], L.sd=0.02)
    ## calculate the logliklihood of each parameter combination
    lik <- lapply(joblist, function(x) {
        allpars <- sort.params(x,fixed.pars,parameter.order);
        estpars <- names(x);
        estpars <- estpars[-which(estpars=='id')];
        trajectory.matching(data=observed.data[[i]], pars=allpars, est=estpars, eval.only=TRUE)
    })
    ## convert results to a dataframe
    results = ldply(lik, function(x) c(x$par,loglik=x$lik))
    ## sort by log-likelihood
    results = results[order(results$loglik),]
    save(results, file=paste0('estimation_results_',i,'_stage_one.rda'))

    ## run subplex on the top 100 parameter sets to hone in more
    tm.results <- vector(mode='list', length=100)
    for (j in 1:100) {
        y <- trajectory.matching(data=observed.data[[i]],
                                 pars=results[j,2:17],
                                 est=estimated.pars,
                                 method='subplex')
        tm.results[[j]] <- y
    }
    tm.results <- ldply(tm.results, function(x) c(x$par, loglik=x$lik))
    tm.results <- tm.results[order(tm.results$loglik),]
    save(tm.results, file=paste0('estimation_results_',i,'_stage_two.rda'))

    ## I am going to take the top 20, create novel variation around each
    ## putative optimal parameter set, and then run subplex again
    tm.vary <- array(NA, dim=c(200,16))
    for (q in 1:20) {
        tm.vary[((q-1)*10+1):(q*10),] <-
            t(sapply(1:10, function(x) {
                n <- names(tm.results[q,1:16])
                w <- which(n%in%estimated.pars)
                p <- as.numeric(tm.results[q,1:16])
                p[w] <- rnorm(length(p[w]),mean=p[w],sd=p[w]/10)
                p
            }))
    }
    colnames(tm.vary) <- colnames(tm.results)[1:16]
    ## Double check to make sure that neither eA nor kappa is greater than 1
    tm.vary[which(tm.vary[,'kappa'] > 1),'kappa'] <- 0.99
    tm.vary[which(tm.vary[,'eA']>1),'eA'] <- 0.99

    tm.vary.results <- vector(mode='list', length=200)
    jvec <- vector()
    for (j in 1:200) {
        jvec <- c(jvec,as.character(j))
        write(jvec, file='progress.txt')
        y <- trajectory.matching(data=observed.data[[i]],
                                     pars=tm.vary[j,],
                                     est=estimated.pars,
                                     method='subplex')
        tm.vary.results[[j]] <- y
        save(tm.vary.results, file=paste0('estimation_results_',i,'_stage_three.rda'))
    }
    tm.vary.results <- ldply(tm.vary.results, function(x) c(x$par, loglik=x$lik, conv=x$conv))
    tm.vary.results <- tm.vary.results[order(tm.vary.results$loglik),]
    save(tm.vary.results, file=paste0('estimation_results_',i,'_stage_three.rda'))

}

