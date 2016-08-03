library(deSolve)
library(magrittr)
library(subplex)
library(pomp)
library(parallel)
library(plyr)
library(tidyr)
library(ggplot2)

## rebuild source for this computer
if (is.loaded("debStochEnv.so")) dyn.unload("debStochEnv.so")
system("rm debStochEnv.so")
system("R CMD SHLIB debStochEnv.c")
dyn.load("debStochEnv.so")

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

optimizer <- function(estpars, fixpars, parorder, transform, obsdata, Np, eval.only=FALSE, type="trajectory_matching", method="subplex") {
    print(estpars)
    estpars <- par_transform(pars=estpars, transform=transform)
    if (any(is.na(estpars)))
        opt <- list(params=estpars,
                    lik=NA)
    else if (eval.only==TRUE) {
        if (type=="trajectory_matching")
            x <- tm_obj(estpars=estpars,
                        data=obsdata,
                        fixpars=fixpars,
                        parorder=parorder,
                        transform=transform)
        else if (type=="particle_filter")
            x <- pf_obj(estpars=estpars,
                        data=obsdata,
                        fixpars=fixpars,
                        parorder=parorder,
                        transform=transform,
                        Np=Np)
        opt <- list(params=estpars,
                    lik=x)
    }
    else {
        if (type=="trajectory_matching") {
            if (method=="subplex")
                x <- subplex(par=estpars,
                             fn=tm_obj,
                             data=obsdata,
                             fixpars=fixpars,
                             parorder=parorder,
                             transform=transform)
            else
                x <- optim(par=estpars,
                           fn=tm_obj,
                           method=method,
                           data=obsdata,
                           fixpars=fixpars,
                           parorder=parorder,
                           transform=transform,
                           control=list(maxit=5000))
        }
        else if (type=="particle_filter") {
            if (method=="subplex")
                x <- subplex(par=estpars,
                             fn=pf_obj,
                             data=obsdata,
                             fixpars=fixpars,
                             parorder=parorder,
                             transform=transform,
                             Np=Np)
            else
                x <- optim(par=estpars,
                           fn=pf_obj,
                           method=method,
                           data=obsdata,
                           fixpars=fixpars,
                           parorder=parorder,
                           transform=transform,
                           Np=Np,
                           control=list(maxit=5000))
        }

        opt <- list(params=par_untransform(x$par,transform),
                    lik=x$value,
                    conv=x$convergence)
    }
    return(opt)
}

## Trajectory matching
tm_obj <- function(estpars, data, fixpars, parorder, transform) {
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

    ## Generic set of food addition events
    eventdat <- data.frame(var="F",
                           time=1:35,
                           value=unname(pars["F0"]),
                           method=rep(c(rep("add",4),"rep"),max(data$age)/5))

    ## Initial condition
    y0 <- c(F=unname(pars["F0"]),
           E=0.00025,
           W=0.00025,
           R=0)

    ## Simulate the system
    try(ode(y0,
        times=0:35,
        func="derivs",
        parms=pars,
        dllname="debStochEnv",
        initfunc="initmod",
        events=list(data=eventdat))) -> out
    if (inherits(out, "try-error"))
        lik <- -Inf
    else if (max(out[,"time"]) < max(data$time))
        lik <- -Inf
    else if (any(is.nan(out)))
        lik <- -Inf
    else if (any(out < 0))
        lik <- -Inf
    ## If no errors, compute the likelihood
    else {
        as.data.frame(out[unique(data$age)+1,]) -> pred

        ## compute the probability of observing the data, given the prediction
        sapply(unique(data$age),
               function(d)
                   c(dnorm(x=data$length[data$age==d],
                           mean=(pred$W[pred$time==d]/pars['xi'])^(1/pars['q']),
                           sd=pars["Lobs"],
                           log=TRUE) %>% sum,
                     dnorm(x=data$eggs[data$age==d],
                           mean=pred$R[pred$time==d],
                           sd=pars["Robs"],
                           log=TRUE) %>% sum
                     ) %>% sum
               ) %>% sum -> lik
    }
    return(-lik)
}


## Implement a particle filter to compute the likelihood of a particular parameter set.
pf_obj <- function(estpars, data, fixpars, parorder, transform, Np=100) {
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

    ## Generic set of food addition events
    eventdat <- data.frame(var="F",
                           time=1:35,
                           value=unname(pars['F0']),
                           method=rep(c(rep("add",4),"rep"),max(data$age)/5))

    ## Generate Np sets of initial conditions to initilize both the filtering distribution and the prediction distribution
    data.frame(T=0,
               F=rnorm(Np, mean=1e6/30, sd=pars['Ferr']),
               E=0.00025,
               W=0.00025,
               R=0) -> x.F -> x.P

    ## timepoints when the likelihood should be evaluated
    times <- c(0, data$age %>% unique)
    tstep <- 1
    lik <- 0
    while (tstep < length(times)) {
        ## For each of the Np particles
        for (i in 1:Np) {
            ## generate the food addition events
            events <- subset(eventdat, time >= times[tstep] & time < times[tstep+1])
            events$value <- rnorm(nrow(events),
                                  mean=events$value,
                                  sd=pars['Ferr'])
            ## obtain a sample of points from the prediction distribution by simulating the model forward
            try(ode(x.F[i,2:5] %>% unlist,
                    times=seq(times[tstep], times[tstep+1]),
                    func="derivs",
                    parms=pars,
                    dllname="debStochEnv",
                    initfunc="initmod",
                    events=list(data=events))) -> out
            if (inherits(out, "try-error") ||
                max(out[,"time"]) < times[tstep+1] ||
                any(is.nan(out)) ||
                any(out < -1e-4)) {
                print(i)
                print(out)
                x.P[i,] <- rep(NA,5)
            }
            else tail(out,1) -> x.P[i,]
        }

        ## determine the weights by computing the probability of observing the data, given the points in the prediction distribution
        sapply((x.P$W/pars['xi'])^(1/pars['q']),
               function(l)
                   dnorm(x=data$length[data$age==times[tstep+1]],
                         mean=l,
                         sd=pars['Lobs'],
                         log=FALSE) %>% sum
               ) +
            sapply(x.P$R,
                   function(r)
                       dnorm(x=data$eggs[data$age==times[tstep+1]],
                             mean=r,
                             sd=pars['Robs'],
                             log=FALSE) %>% sum
                   ) -> weights
        ## set weight to 0 for any particles that had integration errors
        if (any(is.na(weights)))
            weights[is.na(weights)] <- 0

        ## conditional likelihood for this timestep is the mean probability across the points in the prediction distribution
        lik <- lik + log(mean(weights))

        ## use the weights to update the filtering distribution by resampling from the prediction distribution
        w <- cumsum(weights)
        du <- max(w)/length(w)
        u <- runif(1,-du,0)
        p <- vector(mode='numeric', length=length(w))
        for (j in 1:length(w)) {
            i <- 1
            u <- u+du
            while (u > w[i])
                i <- i+1
            p[j] <- i
        }

        x.F <- x.P[p,]
        rownames(x.F) <- as.character(1:length(weights))
        tstep <- tstep+1

    }
    return(-lik)
}
