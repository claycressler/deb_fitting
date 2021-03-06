
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
                x <- try(subplex(par=estpars,
                             fn=tm_obj,
                             data=obsdata,
                             fixpars=fixpars,
                             parorder=parorder,
                             transform=transform))
            else
                x <- try(optim(par=estpars,
                           fn=tm_obj,
                           method=method,
                           data=obsdata,
                           fixpars=fixpars,
                           parorder=parorder,
                           transform=transform,
                           control=list(maxit=5000)))
        }
        else if (type=="particle_filter") {
            if (method=="subplex")
                x <- try(subplex(par=estpars,
                             fn=pf_obj,
                             data=obsdata,
                             fixpars=fixpars,
                             parorder=parorder,
                             transform=transform,
                             Np=Np,
			     control=list(maxit=1000)))
            else
                x <- try(optim(par=estpars,
                           fn=pf_obj,
                           method=method,
                           data=obsdata,
                           fixpars=fixpars,
                           parorder=parorder,
                           transform=transform,
                           Np=Np,
                           control=list(maxit=5000)))
        }
        if (inherits(x, 'try-error'))
            opt <- list(params=par_untransform(estpars,transform),
                        lik=NA,
                        conv=NA)
        else
            opt <- list(params=par_untransform(x$par,transform),
                        lik=x$value,
                        conv=x$convergence)
    }
    return(opt)
}

optimizer2 <- function(estpars, fixpars, parorder, transform, obsdata, Np, eval.only=FALSE, type="trajectory_matching", method="subplex") {
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
                x <- try(subplex(par=estpars,
                             fn=tm_obj,
                             data=obsdata,
                             fixpars=fixpars,
                             parorder=parorder,
                             transform=transform))
            else
                x <- try(optim(par=estpars,
                           fn=tm_obj,
                           method=method,
                           data=obsdata,
                           fixpars=fixpars,
                           parorder=parorder,
                           transform=transform,
                           control=list(maxit=5000)))
        }
        else if (type=="particle_filter") {
            if (method=="subplex")
                x <- try(subplex(par=estpars,
                             fn=pf_obj2,
                             data=obsdata,
                             fixpars=fixpars,
                             parorder=parorder,
                             transform=transform,
                             Np=Np,
			     control=list(maxit=1000)))
            else
                x <- try(optim(par=estpars,
                           fn=pf_obj2,
                           method=method,
                           data=obsdata,
                           fixpars=fixpars,
                           parorder=parorder,
                           transform=transform,
                           Np=Np,
                           control=list(maxit=5000)))
        }
        if (inherits(x, 'try-error'))
            opt <- list(params=par_untransform(estpars,transform),
                        lik=NA,
                        conv=NA)
        else
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

pf_obj2 <- function(estpars, data, fixpars, parorder, transform, Np=100) {
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

    data.frame(F=rnorm(Np, mean=1e6/30, sd=pars['Ferr']),
               E=0.00025,
               W=0.00025,
               R=0) %>% apply(., 1, as.list) %>% lapply(., unlist) -> x0

    weights <- vector(mode='list', length=Np)
    times <- data$age %>% unique
    for (i in 1:Np) {
        eventdat <- data.frame(var="F",
                               time=1:35,
                               value=rnorm(35, mean=unname(pars['F0']), sd=unname(pars['Ferr'])),
                               method=rep(c(rep("add",4),"rep"),max(data$age)/5))
        try(ode(x0[[i]],
                times=0:35,
                func="derivs",
                parms=pars,
                dllname="debStochEnv",
                initfunc="initmod",
                events=list(data=eventdat))) -> out
        if (inherits(out, "try-error") ||
            max(out[,"time"]) < 35 ||
            any(is.nan(out)) ||
            any(out < -1e-4))
            weights[[i]] <- rep(0, length(times))
        else {
            out <- as.data.frame(out)
            ## Calculate the weights for each particle at each observation
            sapply(times,
                   function(t)
                       sapply((out$W[out$time==t]/pars['xi'])^(1/pars['q']) %>% unname,
                              function(w)
                                  dnorm(x=data$length[data$age==t],
                                        mean=w,
                                        sd=pars['Lobs'],
                                        log=FALSE)
                              ) %>% prod *
                                  sapply(out$R[out$time==t] %>% unname,
                                         function(r)
                                             dnorm(x=data$eggs[data$age==t],
                                                   mean=r,
                                                   sd=pars['Robs'],
                                                   log=FALSE) %>% prod
                                         )
                   ) -> weights[[i]]
        }
    }
    ## likelihood of first datapoint
    w <- lapply(weights, function(w) w[1]) %>% unlist
    lik <- log(mean(w[w > 0]))
    for (t in 2:length(times)) {
        ## Begin by resampling the particles based on the weights from the
        ## *previous* timepoint
        v <- cumsum(w)
        du <- max(v)/length(v)
        u <- runif(1,-du,0)
        resamp <- sapply(seq(u+du,u+length(v)*du,du), function(u) min(which(!(u > v))))
        for (i in 1:length(weights)) weights[[i]] <- weights[[resamp[i]]]
        ## Using the weights of these resampled particles, compute the
        ## likelihood for the current timsestep
        w <- lapply(weights, function(w) w[t]) %>% unlist
        lik <- lik + log(mean(w[w > 0]))
    }
    print(-lik)
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
            # obtain a sample of points from the prediction distribution by simulating the model forward
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
                any(out < -1e-4))
                x.P[i,] <- rep(NA,5)
            else tail(out,1) -> x.P[i,]
        }

        ## determine the weights by computing the probability of observing the data, given the points in the prediction distribution
        sapply((x.P$W/pars['xi'])^(1/pars['q']),
               function(l)
                   dnorm(x=data$length[data$age==times[tstep+1]],
                         mean=l,
                         sd=pars['Lobs'],
                         log=FALSE) %>% prod
               ) *
            sapply(x.P$R,
                   function(r)
                       dnorm(x=data$eggs[data$age==times[tstep+1]],
                             mean=r,
                             sd=pars['Robs'],
                             log=FALSE) %>% prod
                   ) -> weights
        ## set weight to 0 for any particles that had integration errors
        if (any(is.na(weights)))
            weights[is.na(weights)] <- 0

        if (any(weights > 0)) {
            ## conditional likelihood for this timestep is the mean probability across the points in the prediction distribution (discarding all zeros)
            lik <- lik + log(mean(weights[weights > 0]))
            ## use the weights to update the filtering distribution by resampling from the prediction distribution
            w <- cumsum(weights)
            du <- max(w)/length(w)
            u <- runif(1,-du,0)
            p <- sapply(seq(u+du,u+length(w)*du,du), function(u) min(which(!(u > w))))
            x.F <- x.P[p,]
            rownames(x.F) <- as.character(1:length(weights))
            tstep <- tstep+1
        }
        else {
	    lik <- -Inf
	    break
        }
    }
    return(-lik)
}

theta <- par_transform(estpars, transform)
sigma <- sqrt(abs(estpars)/10)

it_filter <- function(theta, sigma, data, fixpars, parorder, transform, cf=0.98, vf=2, B=1000, N=100) {
    ## set up storage for the parameters through time
    theta.k <- array(NA, dim=c(length(estpars), N+1))
    ## initialize the parameters
    theta.k[,1] <- theta
    ## set up storage for the random parameters
    theta.i <- array(NA, dim=c(length(estpars), B))
    rownames(theta.k) <- rownames(theta.i) <- names(theta)
    ## set up storage for the state variables
    states.i <- array(NA, dim=c(4, B))
    rownames(states.i) <- c("F","E","W","R")
    ## set up storage for the weights
    weights <- vector(mode='numeric', length=B)

    ## Generic set of food addition events
    eventdat <- data.frame(var="F",
                           time=1:35,
                           value=unname(pars['F0']),
                           method=rep(c(rep("add",4),"rep"),max(data$age)/5))

    for (n in 2:(N+1)) { ## N is the number of iterations
        ## draw an initial set of B particles and states
        for (b in 1:B) {
            theta.i[,b] <- rnorm(nrow(theta.i),mean=theta.k[,n-1],sd=vf*cf*sigma)
            states.i[,b] <- c(1e6/30, 0.00025, 0.00025, 0)
        }
        ## do particle filtering on the random-parameters model
        times <- c(0, unique(data$age)); tstep <- 1
        while (tstep < length(times)) {
            ## for each of the particles
            for (b in 1:B) {
                ## Put the parameters back on the natural scale
                estpars <- par_untransform(theta.i[,b], transform)
                ## compute Imax and g based on the estimate of fh
                fixpars["Imax"] <- calc_Imax(unname(estpars["fh"]))
                fixpars["g"] <- calc_g(unname(estpars["fh"]))
                ## combine the parameters to be estimated and the fixed parameters
                ## into a single vector, with an order specified by parorder
                pars <- c(estpars, fixpars)
                pars[match(parorder, names(pars))] -> pars
                if(any(names(pars)!=parorder)) stop("parameter are in the wrong order")
                ## Randomize the food additions
                events <- subset(eventdat, time >= times[tstep] & time < times[tstep+1])
                events$value <- rnorm(nrow(events),
                                      mean=events$value,
                                      sd=pars['Ferr'])
                ## simulate the system one timestep forward
                try(ode(states.i[,b],
                        times=seq(times[tstep], times[tstep+1]),
                        func="derivs",
                        parms=pars,
                        dllname="debStochEnv",
                        initfunc="initmod",
                        events=list(data=events))) -> out
                ## if there were no errors, compute the weight for this particle
                if (inherits(out, "try-error") || max(out[,"time"]) < times[tstep+1] ||
                    any(is.nan(out)) || any(out < -1e-4))
                    weights[b] <- 0
                else
                    sum(dnorm(x=data$length[data$age==times[tstep+1]],
                          mean=tail(out[,"W"],1),
                          sd=pars['Lobs'],
                          log=TRUE) +
                              dnorm(x=data$eggs[data$age==times[tstep+1]],
                                    mean=tail(out[,"R"],1),
                                    sd=pars['Robs'],
                                    log=TRUE)) -> weights[b]
            }



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


