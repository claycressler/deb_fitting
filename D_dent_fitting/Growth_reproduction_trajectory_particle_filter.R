source("Growth_reproduction_trajectory_fitting_stochastic_functions.R")

pars <- c(Imax=22500, fh=10000, g=1.45, rho=0.1, eps=44.5e-9, V=30, F0=1000000/30, xi=2.62e-3, q=2.4, K=0.3, km=0.15, ER=1.51e-3, v=10, Lobs=0.1, Robs=1, Ferr=10000)
y0 <- c(F=1000000/30, E=0.00025, W=0.00025, R=0)

## days when observations will take place
times <- c(5,10,12,15,18,25,30,35)
## number of reps per observation day
inds <- 1:12

set.seed(101)
out <- vector(mode='list', length=length(times)*length(inds))
for (i in 1:(length(times)*length(inds))) {
    ## feeding schedule - amount of food added each time is stochastic
    eventdat <- data.frame(var="F",
                           time=1:35,
                           value=rnorm(35,
                               mean=unname(pars["F0"]),
                               sd=unname(pars["Ferr"])),
                           method=rep(c(rep("add",4),"rep"),max(times)/5))
    ode(y0, times=0:35, func="derivs", parms=pars, dllname="debStochEnv", initfunc="initmod", events=list(data=eventdat)) %>% as.data.frame -> out[[i]]
}

## Simulate observations
set.seed(101)
data = expand.grid(ind=inds, age=times)
data$length <- rnorm(nrow(data),
                     mean=(sapply(1:nrow(data), function(i) out[[i]][data[i,'age']+1,4])/pars['xi'])^(1/pars['q']),
                     sd=pars['Lobs'])
data$eggs <- rnorm(nrow(data),
                   mean=sapply(1:nrow(data), function(i) out[[i]][data[i,'age']+1,5]),
                   sd=pars['Robs'])

fixpars <- pars[c("Imax","g","eps","V","F0","xi","q","ER","v")]
estpars <- pars[c("fh","rho","K","km","Lobs","Robs","Ferr")]
transform <- c("log", rep("logit",2), rep("log",4))
parorder <- c("Imax","fh","g","rho","eps","V","F0","xi","q","K","km","ER","v","Lobs","Robs","Ferr")

pf_obj(par_transform(estpars,transform), data, fixpars, parorder, transform, 100)
optimizer(estpars, fixpars, parorder, transform, obsdata=data, Np=100, eval.only=TRUE, type="particle_filter")

tm_obj(par_transform(estpars,transform), data, fixpars, parorder, transform)
optimizer(estpars, fixpars, parorder, transform, obsdata=data, eval.only=TRUE, type="trajectory_matching")

######################################################################
######################################################################
######################################################################

## Generate 25 different parameter sets and use each to generate simulated data
datasets <- vector(mode='list', length=25)
set.seed(555)
for (d in 1:25) {
    pars <- c(Imax=22500,
              fh=rnorm(1,mean=10000,sd=1000),
              g=1.45,
              rho=exp(rnorm(1,mean=log(0.1),sd=1)),
              eps=44.5e-9,
              V=30,
              F0=1000000/30,
              xi=2.62e-3,
              q=2.4,
              K=rnorm(1,mean=0.3,sd=0.1),
              km=rnorm(1,mean=0.15,sd=0.05),
              ER=rnorm(1,mean=1.5e-3,sd=2.5e-4),
              v=10,
              Lobs=rnorm(1, mean=0.2, sd=0.04),
              Robs=rnorm(1, mean=3, sd=0.5),
              Ferr=rnorm(1, mean=10000, sd=5000))

    y0 <- c(F=1000000/30, E=0.00025, W=0.00025, R=0)

    ## days when observations will take place
    times <- c(5,10,12,15,18,25,30,35)
    ## number of reps per observation day
    inds <- 1:12

    out <- vector(mode='list', length=length(times)*length(inds))
    for (i in 1:(length(times)*length(inds))) {
        ## feeding schedule - amount of food added each time is stochastic
        eventdat <- data.frame(var="F",
                               time=1:35,
                               value=rnorm(35,
                                   mean=unname(pars["F0"]),
                                   sd=unname(pars["Ferr"])),
                               method=rep(c(rep("add",4),"rep"),max(times)/5))
        ode(y0, times=0:35, func="derivs", parms=pars, dllname="debStochEnv", initfunc="initmod", events=list(data=eventdat)) %>% as.data.frame -> out[[i]]
    }

    ## Simulate observations
    data = expand.grid(ind=inds, age=times)
    data$length <- rnorm(nrow(data),
                         mean=(sapply(1:nrow(data), function(i) out[[i]][data[i,'age']+1,4])/pars['xi'])^(1/pars['q']),
                         sd=pars['Lobs'])
    data$eggs <- rnorm(nrow(data),
                       mean=sapply(1:nrow(data), function(i) out[[i]][data[i,'age']+1,5]),
                       sd=pars['Robs'])
    datasets[[d]] <- list(params=pars, data=data)
}

saveRDS(datasets, file="env_stoch_datasets.RDS")

pars <- c(Imax=22500, fh=10000, g=1.45, rho=0.1, eps=44.5e-9, V=30, F0=1000000/30, xi=2.62e-3, q=2.4, K=0.3, km=0.15, ER=1.51e-3, v=10, Lobs=0.1, Robs=1, Ferr=10000)
fixpars <- pars[c("Imax","g","eps","V","F0","xi","q","ER","v")]
estpars <- pars[c("fh","rho","K","km","Lobs","Robs","Ferr")]
transform <- c("log", rep("logit",2), rep("log",4))
parorder <- c("Imax","fh","g","rho","eps","V","F0","xi","q","K","km","ER","v","Lobs","Robs","Ferr")


## For this first set of attempts, do not attempt to estimate ER, but allow it to be fixed at the correct value.
tm_ests <- vector(mode='list', length=25)
for (d in 2:25) {
    print(d)
    data <- datasets[[d]]$data

    ## Begin the funnel of optimization with simple trajectory matching,
    ## then use the particle filter to hone the parameter estimates.
    box <- cbind(lower=c(fh=2000, rho=0, K=0, km=0.001, Lobs=0.001, Robs=0.01, Ferr=1000),
                 upper=c(fh=20000, rho=1, K=1, km=1, Lobs=2, Robs=10, Ferr=50000))
    sobolDesign(lower=box[,'lower'],
                upper=box[,'upper'],
                nseq=20000) %>%
                    apply(., 1, as.list) %>%
                        lapply(., unlist) -> guesses

    fixpars <- c(Imax=22500, g=1.45, eps=44.5e-9, V=30, F0=1000000/30, xi=2.62e-3, q=2.4, ER=unname(datasets[[d]]$params['ER']), v=100)
    transform <- c("log", rep("logit",2), rep("log",4))
    parorder <- c("Imax","fh","g","rho","eps","V","F0","xi","q","K","km","ER","v","Lobs","Robs","Ferr")

    mclapply(guesses,
             optimizer,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=data,
             eval.only=TRUE,
             type='trajectory_matching',
             mc.cores=15) %>%
                 lapply(., function(x) x$lik) %>%
                     unlist -> guess_lik
    guesses[order(guess_lik)[1:300]] -> refine
    mclapply(refine,
             optimizer,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=data,
             eval.only=FALSE,
             type='trajectory_matching',
             method='Nelder-Mead',
             mc.cores=15) -> refine_lik
    refine_lik %>%
        lapply(., unlist) %>%
            unlist %>%
                matrix(., ncol=(nrow(box)+2), byrow=TRUE, dimnames=list(1:length(refine_lik), c(rownames(box),"lik","conv"))) %>%
                    as.data.frame -> refine_pars
    refine_pars[order(refine_pars$lik),] -> refine_pars
    tm_ests[[d]] <- refine_pars
    saveRDS(tm_ests, file="Trajectory_matching_7-27.RDS")
}

saveRDS(datasets, file="env_stoch_datasets.RDS")

## For this first set of attempts, do not attempt to estimate ER, but allow it to be fixed at the correct value.
source("Growth_reproduction_trajectory_fitting_stochastic_functions.R")
pf_ests <- vector(mode='list', length=25)
tm_ests <- readRDS("Trajectory_matching_7-27.RDS")
datasets <- readRDS("env_stoch_datasets.RDS")

pars <- c(Imax=22500, fh=10000, g=1.45, rho=0.1, eps=44.5e-9, V=30, F0=1000000/30, xi=2.62e-3, q=2.4, K=0.3, km=0.15, ER=1.51e-3, v=10, Lobs=0.1, Robs=1, Ferr=10000)
fixpars <- pars[c("Imax","g","eps","V","F0","xi","q","ER","v")]
estpars <- pars[c("fh","rho","K","km","Lobs","Robs","Ferr")]
transform <- c("log", rep("logit",2), rep("log",4))
parorder <- c("Imax","fh","g","rho","eps","V","F0","xi","q","K","km","ER","v","Lobs","Robs","Ferr")

for (d in 5:25) {
    print(d)
    tic <- Sys.time()
    print(tic)
    data <- datasets[[d]]$data
    estpars <- tm_ests[[d]]
    estpars$Ferr <- 1e4
    ## pull out the first 15 distinct parameter sets and use those in the particle filter (because many of the parameter sets are very similar to one another)
    inds <- 1; i <- 1
    while (length(inds) < 15) {
        ## check to see if this fh value is distinct from all others currently in the set
        if(!any(signif(estpars$rho[i],2)==signif(estpars$rho[inds],2)))
            inds <- c(inds, i)
        i <- i+1
    }
    lapply(apply(estpars[inds,1:7], 1, as.list), unlist) -> refine
    mclapply(refine,
             optimizer,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=data,
             eval.only=FALSE,
             type='particle_filter',
             method='subplex',
             Np=100,
             mc.cores=15) -> refine_lik
    refine_lik %>%
        lapply(., unlist) %>%
            unlist %>%
                matrix(., ncol=ncol(estpars)+1, byrow=TRUE, dimnames=list(1:length(refine_lik), c(colnames(estpars),'error'))) %>%
                    as.data.frame -> refine_pars
    refine_pars[order(refine_pars$lik),] -> refine_pars
    pf_ests[[d]] <- refine_pars
    toc <- Sys.time()
    print(tic); print(toc-tic)
    saveRDS(pf_ests, file="Particle_filter_7-30.RDS")
}











## for this set, allow it to be estimated along with everything else

tm_ests <- vector(mode='list', length=25)
for (d in 1:25) {
    print(d)
    data <- datasets[[d]]$data

    ## Begin the funnel of optimization with simple trajectory matching,
    ## then use the particle filter to hone the parameter estimates.
    box <- cbind(lower=c(fh=2000, rho=0, K=0, km=0.001, ER=1e-4, Lobs=0.001, Robs=0.01, Ferr=1000),
                 upper=c(fh=20000, rho=1, K=1, km=1, ER=1e-2, Lobs=2, Robs=10, Ferr=50000))
    sobolDesign(lower=box[,'lower'],
                upper=box[,'upper'],
                nseq=20000) %>%
                    apply(., 1, as.list) %>%
                        lapply(., unlist) -> guesses

    fixpars <- c(Imax=22500, g=1.45, eps=44.5e-9, V=30, F0=1000000/30, xi=2.62e-3, q=2.4, v=100)
    transform <- c("log", rep("logit",2), rep("log",5))
    parorder <- c("Imax","fh","g","rho","eps","V","F0","xi","q","K","km","ER","v","Lobs","Robs","Ferr")

    mclapply(guesses,
             optimizer,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=data,
             eval.only=TRUE,
             type='trajectory_matching',
             mc.cores=15) %>%
                 lapply(., function(x) x$lik) %>%
                     unlist -> guess_lik
    guesses[order(guess_lik)[1:300]] -> refine
    mclapply(refine,
             optimizer,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=data,
             eval.only=FALSE,
             type='trajectory_matching',
             method='Nelder-Mead',
             mc.cores=15) -> refine_lik
    refine_lik %>%
        lapply(., unlist) %>%
            unlist %>%
                matrix(., ncol=(nrow(box)+2), byrow=TRUE, dimnames=list(1:length(refine_lik), c(rownames(box),"lik","conv"))) %>%
                    as.data.frame -> refine_pars
    refine_pars[order(refine_pars$lik),] -> refine_pars
    tm_ests[[d]] <- refine_pars
    saveRDS(tm_ests, file="Trajectory_matching_7-27.RDS")
}




