library(deSolve)
library(magrittr)
library(subplex)
library(pomp)
library(parallel)
library(plyr)

## rebuild source for this computer
system("rm deb.so")
system("R CMD SHLIB deb.c")
dyn.load("deb.so")

## Transform the parameters.  The optimization routines work best if
## the parameter values are unconstrained. Most of the parameters vary
## between 0 and Inf; these can be made unconstrained by simple
## log-transforming. The parameter K, on the other hand, can only vary
## between 0 and 1. To put a parameter that varies between 0 and 1 on
## an unconstrained scale, we use the logit transform.
par_transform <- function(pars, transform) {
    for (i in 1:length(pars)) {
        if (transform[i]=="log")
            pars[i] <- log(pars[i])
        if (transform[i]=="logit")
            pars[i] <- log(pars[i]/(1-pars[i]))
    }
    return(pars)
}
## Untransform the parameters
par_untransform <- function(pars, transform) {
    for (i in 1:length(pars)) {
        if (transform[i]=="log")
            pars[i] <- exp(pars[i])
        if (transform[i]=="logit")
            pars[i] <- exp(pars[i])/(1+exp(pars[i]))
    }
    return(pars)
}
traj_match <- function(estpars, fixpars, parorder, transform, obsdata, eval.only=FALSE, method="subplex") {
    estpars <- par_transform(pars=estpars, transform=transform)
    if (any(is.na(estpars)))
        opt <- list(params=estpars,
                    lik=NA)
    else if (eval.only==TRUE) {
        x <- obj(estpars=estpars,
                     data=obsdata,
                     fixpars=fixpars,
                     parorder=parorder,
                     transform=transform)
        opt <- list(params=estpars,
                    lik=x)
    }
    else {
        if (method=="subplex")
            x <- subplex(par=estpars,
                         fn=obj,
                         data=obsdata,
                         fixpars=fixpars,
                         parorder=parorder,
                         transform=transform)
        else
            x <- optim(par=estpars,
                       fn=obj,
                       method=method,
                       data=obsdata,
                       fixpars=fixpars,
                       parorder=parorder,
                       transform=transform,
                       control=list(maxit=5000))
        opt <- list(params=par_untransform(x$par,transform),
                    lik=x$value,
                    conv=x$convergence)
    }
    return(opt)
}
## Objective function to minimize for the structured model
obj <- function(estpars, data, fixpars, parorder, transform) {
    ## We will give the model the true initial conditions
    y0 <- c(E=0.0002, W=0.0005, R=0)
    ## Put the estimated parameters back on the natural scale
    estpars <- par_untransform(estpars, transform)
    ## combine the parameters to be estimated and the fixed parameters
    ## into a single vector, with an order specified by parorder
    pars <- c(estpars, fixpars)
    pars[match(parorder, names(pars))] -> pars
    if(any(names(pars)!=parorder)) stop("parameter are in the wrong order")
    ## simulate the model at these parameters
    try(ode(y0,
            times=seq(0, max(data$time)),
            func="derivs",
            parms=pars,
            dllname="deb",
            initfunc="initmod"
            )) -> out
    ## if this parameter set produces a simulation error, if it fails
    ## to completely run, or if any of the variables take negative
    ## values, return -Inf for the likelihood.
    if (inherits(out, "try-error"))
        lik <- -Inf
    else if (max(out[,"time"]) < max(data$time))
        lik <- -Inf
    else if (any(is.nan(out)))
        lik <- -Inf
    else if (any(out < 0))
        lik <- -Inf
    else {
        out <- as.data.frame(out)
        mutate(out,
               L=((W+E)/pars["xi"])^(1/pars["q"]),
               R=R-R[which((E+W) < 5.9e-3) %>% max]) -> out
        out$R[out$R < 0] <- 0
        ## calculate the log-likelihood of observing these data
        ## assuming only observation error
        (dnorm(data$length,
              sapply(data$times,
                     function(t) out$L[out$time==t]),
              unname(pars["Lobs"]),
              log=TRUE) +
        dpois(data$eggs,
              sapply(data$times,
                     function(t) (out$R-out$R[which(out$E+out$W < 5.9e-3) %>% max])[t+1]),
              log=TRUE)) %>%
            sum(., na.rm=TRUE) -> lik
    }
    ## return the negative log-likelihood (because optimization methods minimize)
    return(-lik)
}

pars <- c(rho=0.7, eps=40e-9, Imax=0.0126, g=1.37, F=1e5/0.01, xi=2.62e-3, q=2.4, K=0.6, km=0.33, ER=1.51e-3, v=10, Lobs=0.1)
ode(y=c(E=0.0002, W=0.0005, R=0),
    times=seq(0, 35),
    func="derivs",
    parms=pars,
    dllname="deb",
    initfunc="initmod"
    ) %>% as.data.frame -> out
## assume that length is observed assuming a normal distribution
## total; observed eggs is Poisson-distributed, but subtract off all
## "reproduction" done prior to reaching sexual maturity when total
## mass is greater than 5.9e-3 mg (Nisbet et al. 2004).
days <- c(5,10,12,15,18,25,30,35)
set.seed(1239478)
data.frame(times=rep(days, each=12),
           length=sapply(with(out, ((E[days+1]+W[days+1])/pars["xi"])^(1/pars["q"])),
               function(x)
                   rnorm(12, mean=x, sd=pars["Lobs"])
                         ) %>% as.numeric,
           eggs=sapply(with(out, R-R[which((E+W) < 5.9e-3) %>% max])[days+1],
               function(x)
                   rpois(12, lambda=x)
                       ) %>% as.numeric
           ) -> data

parorder=c("rho","eps","Imax","g","F","xi","q","K","km","ER","v","Lobs")
estpars=c(K=0.6, km=0.33, ER=1.51e-3, v=10, Lobs=0.1)
fixpars=c(rho=0.7, eps=40e-9, Imax=0.0126, g=1.37, F=1e5/0.01, xi=2.62e-3, q=2.4)
transform=c("logit",rep("log",4))
## TESTING
all(par_untransform(par_transform(estpars, transform), transform)-estpars < 1e-10)

ode(y=c(E=0.0002, W=0.0005, R=0),
    times=seq(0, 35),
    func="derivs",
    parms=c(fixpars,estpars),
    dllname="deb",
    initfunc="initmod"
    ) %>% as.data.frame -> out
mutate(out, L=((W+E)/pars["xi"])^(1/pars["q"])) -> out
## calculate the log-likelihood of observing these data
## assuming only observation error
(dnorm(data$length,
       sapply(data$times,
              function(t) out$L[out$time==t]),
       unname(pars["Lobs"]),
       log=TRUE) +
     dpois(data$eggs,
           sapply(data$times,
                  function(t) (out$R-out$R[which(out$E+out$W < 5.9e-3) %>% max])[t+1]),
           log=TRUE)) %>%
    sum(., na.rm=TRUE) -> lik1

obj(par_transform(estpars, transform), data, fixpars, parorder, transform) -> lik2

traj_match(estpars, fixpars, parorder, transform, data, eval.only=TRUE)$lik -> lik3

## GENERATE RANDOM PARAMETER SETS AND OBSERVED DATASETS
parorder=c("rho","eps","Imax","g","F","xi","q","K","km","ER","v","Lobs")
estpars0=c(K=0.6, km=0.33, ER=1.51e-3, v=10, Lobs=0.1)
fixpars=c(rho=0.7, eps=40e-9, Imax=0.0126, g=1.37, F=1e5/0.01, xi=2.62e-3, q=2.4)
transform=c("logit",rep("log",4))

days <- c(5,10,12,15,18,25,30,35)
datasets <- vector(mode='list', length=20)
set.seed(123407)
for (i in 1:20) {
    data.frame(times=rep(days, each=12),
               length=rep(0,96),
               eggs=rep(0,96)) -> data
    while ( (min(data$length) < 0.5) | (max(data$length) > 4) | (max(data[,"eggs"]) < 20) | any(is.na(data)) ) {
        ## GENERATE NOVEL PARAMETERS
        rnorm(length(estpars0), mean=estpars0, sd=estpars0/2) -> p
        names(p) <- names(estpars0)
        while (p["K"] > 0.9 | p["K"] < 0.1 | any(p < 0)) {
            rnorm(length(estpars0), mean=estpars0, sd=estpars0/2) -> p
            names(p) <- names(estpars0)
        }
        estpars <- p
        ## SIMULATE
        ode(y=c(E=0.0002, W=0.0005, R=0),
            times=seq(0, 35),
            func="derivs",
            parms=c(fixpars, estpars),
            dllname="deb",
            initfunc="initmod"
            ) %>% as.data.frame -> out
        mutate(out, R=R-R[which((E+W) < 5.9e-3) %>% max]) -> out
        out$R[out$R < 0] <- 0
        ## SAMPLE TO CREATE A DATASET
        data.frame(times=rep(days, each=12),
                   length=sapply(with(out, ((E[days+1]+W[days+1])/fixpars["xi"])^(1/fixpars["q"])),
                       function(x)
                           rnorm(12, mean=x, sd=estpars["Lobs"])
                                 ) %>% as.numeric,
                   eggs=sapply(out$R[days+1],
                       function(x)
                           rpois(12, lambda=x)
                               ) %>% as.numeric
                   ) -> data
    }
    datasets[[i]] <- list(params=c(fixpars, estpars), data=data)
}


run <- FALSE
if (run) {
    ## just do the first parameter set (see notes for explanation why)
    est_params <- vector(mode='list', length=20)
    data <- datasets[[1]]$data

    ## Generate a large number of different initial parameter guesses
    box <- cbind(lower=c(K=0, km=0.001, ER=0.00001, v=0.1, Lobs=0.0001),
                 upper=c(K=1, km=10, ER=1, v=1000, Lobs=2))
    sobolDesign(lower=box[,'lower'],
                upper=box[,'upper'],
<<<<<<< HEAD
                nseq=100000) %>%
                    apply(., 1, as.list) %>%
                        lapply(., unlist) -> guesses
=======
                nseq=500000) %>%
        apply(., 1, as.list) %>%
            lapply(., unlist) -> guesses
>>>>>>> a0501118b86fb1d6849ef13b54a0a59a2deb2e7a
    mclapply(guesses,
             traj_match,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=data,
             eval.only=TRUE,
<<<<<<< HEAD
             mc.cores=5) %>%
                 lapply(., function(x) x$lik) %>%
                     unlist -> guess_lik
=======
             mc.cores=10) %>%
        lapply(., function(x) x$lik) %>%
            unlist -> guess_lik
>>>>>>> a0501118b86fb1d6849ef13b54a0a59a2deb2e7a

    guesses[order(guess_lik)[1:1000]] -> refine
    mclapply(refine,
             traj_match,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=data,
             eval.only=FALSE,
             method="Nelder-Mead",
             mc.cores=5) -> refine_lik
    refine_lik %>%
        lapply(., function(l) c(l$params, l$lik)) %>%
            unlist %>%
                matrix(., ncol=6, byrow=TRUE) %>%
                    as.data.frame -> refine_pars
    colnames(refine_pars) <- c(rownames(box), "lik")
    est_params[[i]] <-refine_pars
<<<<<<< HEAD
    saveRDS(est_params, file="~/Dropbox/Growth_reproduction_trajectory_fitting_2.RDS")
}
=======
    saveRDS(est_params, file="~/Dropbox/Growth_reproduction_trajectory_fitting.RDS")
>>>>>>> a0501118b86fb1d6849ef13b54a0a59a2deb2e7a

    ## For this run, using the data only from the first dataset, what happens if you fix the estimate of v at various values, and then try to estimate the other parameters?
    data <- datasets[[1]]$data
    parorder=c("rho","eps","Imax","g","F","xi","q","K","km","ER","v","Lobs")
    transform=c("logit",rep("log",3))
    v_vals <- seq(20,2000,20)
    est_params <- vector(mode='list', length=length(v_vals))
    box <- cbind(lower=c(K=0, km=0.001, ER=1e-6, Lobs=0.01),
                 upper=c(K=1, km=1, ER=1e-3, Lobs=1))
    sobolDesign(lower=box[,'lower'],
                upper=box[,'upper'],
                nseq=50000) %>%
        apply(., 1, as.list) %>%
            lapply(., unlist) -> guesses
    for (i in 1:length(v_vals)) {
        v <- v_vals[i]
        print("Growth_reproduction_trajectory_fitting_profile_lik_v")
        print(v)
        fixpars=c(rho=0.7, eps=40e-9, Imax=0.0126, g=1.37, F=1e5/0.01, xi=2.62e-3, q=2.4, v=v)
        t1 <- Sys.time()
        mclapply(guesses,
                 traj_match,
                 fixpars=fixpars,
                 parorder=parorder,
                 transform=transform,
                 obsdata=data,
                 eval.only=TRUE,
                 mc.cores=10) %>%
            lapply(., function(x) x$lik) %>%
                unlist -> guess_lik
        t2 <- Sys.time()
        print(paste0("Initial guesses took", t2-t1))
        guesses[order(guess_lik)[1:100]] -> refine
        t1 <- Sys.time()
        mclapply(refine,
                 traj_match,
                 fixpars=fixpars,
                 parorder=parorder,
                 transform=transform,
                 obsdata=data,
                 eval.only=FALSE,
                 method="Nelder-Mead",
                 mc.cores=10) -> refine_lik
        t2 <- Sys.time()
        print(paste0("Nelder-Mead took", t2-t1))
        refine_lik %>%
            unlist %>%
                matrix(., ncol=6, byrow=TRUE) %>%
                    as.data.frame -> refine_pars
        colnames(refine_pars) <- c(rownames(box), "lik", "conv")
        est_params[[i]] <-refine_pars
        saveRDS(est_params, file="Growth_reproduction_trajectory_fitting_profile_lik_v.RDS")
    }
}

## zoom in on the lower values of v
    data <- datasets[[1]]$data
    parorder=c("rho","eps","Imax","g","F","xi","q","K","km","ER","v","Lobs")
    transform=c("logit",rep("log",3))
    v_vals <- seq(1,20)
    est_params <- vector(mode='list', length=length(v_vals))
    box <- cbind(lower=c(K=0, km=0.001, ER=1e-6, Lobs=0.01),
                 upper=c(K=1, km=1, ER=1e-3, Lobs=1))
    sobolDesign(lower=box[,'lower'],
                upper=box[,'upper'],
                nseq=50000) %>%
        apply(., 1, as.list) %>%
            lapply(., unlist) -> guesses
    for (i in 1:length(v_vals)) {
        v <- v_vals[i]
        print("Growth_reproduction_trajectory_fitting_profile_lik_v")
        print(v)
        fixpars=c(rho=0.7, eps=40e-9, Imax=0.0126, g=1.37, F=1e5/0.01, xi=2.62e-3, q=2.4, v=v)
        t1 <- Sys.time()
        mclapply(guesses,
                 traj_match,
                 fixpars=fixpars,
                 parorder=parorder,
                 transform=transform,
                 obsdata=data,
                 eval.only=TRUE,
                 mc.cores=10) %>%
            lapply(., function(x) x$lik) %>%
                unlist -> guess_lik
        t2 <- Sys.time()
        print(paste0("Initial guesses took", t2-t1))
        guesses[order(guess_lik)[1:100]] -> refine
        t1 <- Sys.time()
        mclapply(refine,
                 traj_match,
                 fixpars=fixpars,
                 parorder=parorder,
                 transform=transform,
                 obsdata=data,
                 eval.only=FALSE,
                 method="Nelder-Mead",
                 mc.cores=10) -> refine_lik
        t2 <- Sys.time()
        print(paste0("Nelder-Mead took", t2-t1))
        refine_lik %>%
            unlist %>%
                matrix(., ncol=6, byrow=TRUE) %>%
                    as.data.frame -> refine_pars
        colnames(refine_pars) <- c(rownames(box), "lik", "conv")
        est_params[[i]] <-refine_pars
        saveRDS(est_params, file="Growth_reproduction_trajectory_fitting_profile_lik_v_refine.RDS")
    }

