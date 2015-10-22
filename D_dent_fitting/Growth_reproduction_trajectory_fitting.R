library(deSolve)
library(magrittr)
library(subplex)
library(pomp)
library(parallel)
library(plyr)

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
    print(lik)
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
        rnorm(length(estpars), mean=estpars0, sd=estpars0/2) -> p
        names(p) <- names(estpars)
        while (p["K"] > 0.9 | p["K"] < 0.1 | any(p < 0)) {
            rnorm(length(estpars), mean=estpars, sd=estpars/2) -> p
            names(p) <- names(estpars)
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
                   length=sapply(with(out, ((E[days+1]+W[days+1])/pars["xi"])^(1/pars["q"])),
                       function(x)
                           rnorm(12, mean=x, sd=pars["Lobs"])
                                 ) %>% as.numeric,
                   eggs=sapply(out$R[days+1],
                       function(x)
                           rpois(12, lambda=x)
                               ) %>% as.numeric
                   ) -> data
    }
    datasets[[i]] <- list(params=c(fixpars, estpars), data=data)
}

est_params <- vector(mode='list', length=20)
for (i in 1:20) {
    data <- datasets[[i]]$data

    ## Generate a large number of different initial parameter guesses
    box <- cbind(lower=c(K=0, km=0.001, ER=0.00001, v=0.1, Lobs=0.0001),
                 upper=c(K=1, km=10, ER=1, v=1000, Lobs=2))
    sobolDesign(lower=box[,'lower'],
                upper=box[,'upper'],
                nseq=500000) %>%
                    apply(., 1, as.list) %>%
                        lapply(., unlist) -> guesses
    mclapply(guesses,
             traj_match,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=data,
             eval.only=TRUE,
             mc.cores=5) %>%
                 lapply(., function(x) x$lik) %>%
                     unlist -> guess_lik

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
    saveRDS(est_params, file="~/Dropbox/Growth_reproduction_trajectory_fitting.RDS")
}



