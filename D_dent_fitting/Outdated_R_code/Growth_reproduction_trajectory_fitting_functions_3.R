library(deSolve)
library(magrittr)
library(subplex)
library(pomp)
library(parallel)
library(plyr)
library(tidyr)
library(ggplot2)

## rebuild source for this computer
if (is.loaded("deb2.so")) dyn.unload("deb2.so")
system("rm deb2.so")
system("R CMD SHLIB deb2.c")
dyn.load("deb2.so")

calc_Imax <- function(fh) 10616.500976 + 1.129421*fh
calc_g <- function(fh) 1.38413 + 9.839020e-6*fh - 2.738144e-10*fh^2 + 2.648817e-15*fh^3

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

traj_match <- function(estpars, fixpars, parorder, transform, obsdata, events, eval.only=FALSE, method="subplex") {
    estpars <- par_transform(pars=estpars, transform=transform)
    if (any(is.na(estpars)))
        opt <- list(params=estpars,
                    lik=NA)
    else if (eval.only==TRUE) {
        x <- obj(estpars=estpars,
                 data=obsdata,
                 fixpars=fixpars,
                 parorder=parorder,
                 transform=transform,
                 events=events)
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
                         transform=transform,
                         events=events)
        else
            x <- optim(par=estpars,
                       fn=obj,
                       method=method,
                       data=obsdata,
                       fixpars=fixpars,
                       parorder=parorder,
                       transform=transform,
                       events=events,
                       control=list(maxit=5000))
        opt <- list(params=par_untransform(x$par,transform),
                    lik=x$value,
                    conv=x$convergence)
    }
    return(opt)
}
## Objective function to minimize for the structured model


obj <- function(estpars, data, fixpars, parorder, transform, events) {
    ## We will give the model the true initial conditions
    y0 <- c(F=1e6/30, E=0.00025, W=0.00025, R=0)
    ## Put the estimated parameters back on the natural scale
    estpars <- par_untransform(estpars, transform)
    ## compute Imax and g based on the estimate of fh
    fixpars["Imax"] <- calc_Imax(unname(estpars["fh"]))
    fixpars["g"] <- calc_g(unname(estpars["fh"]))
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
            dllname="deb2",
            initfunc="initmod",
            events=list(data=events)
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
               R=R-R[which((E+W) < 0.005) %>% max]) -> out
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
                          function(t) out$R[out$time==t]),
                   log=TRUE)) %>%
            sum(., na.rm=TRUE) -> lik
    }
    ## return the negative log-likelihood (because optimization methods minimize)
    return(-lik)
}
