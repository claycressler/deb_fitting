library(deSolve)
library(magrittr)
library(subplex)
library(pomp)
library(parallel)
library(dplyr)

## Transform the parameters.
## The optimization routines work best if the parameter values are
## unconstrained. Most of the parameters vary between 0 and Inf; these
## can be made unconstrained by simple log-transforming. The parameter
## t_rep, on the other hand, can only vary between 0 and the max time
## in the fitted dataset. Thus, we let t_rep take only values between
## 0 and 1, and then multiply the value of t_rep by the max time to
## get the value of t_rep for the simulations. To put a parameter that
## varies between 0 and 1 on an unconstrained scale, we use the logit
## transform.
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
## Trajectory matching to estimate the parameters of the dynamical model.
## estpars is a named numeric vector giving initial guesses for the
##      parameters whose values will be inferred using trajecory matching
## fixpars is a named numeric vector giving the values of parameters
##      whose values are not estimated
## parorder is a character vector giving the order of parameters
##     expected by the C function that simulates the dynamical model
## transform is a character vector giving the transforms for each of the
##     estimated parameters
## obsdata is a dataframe containing the observed data (times and spore counts)
## eval.only is a logical variable that determines whether the
##     likelihood is to maximized or simply calculated for a set of
##     parameters
## method is the name of the optimizer to use
traj_match <- function(estpars, fixpars, parorder, transform, obsdata, eval.only=FALSE, method="subplex") {
    estpars <- par_transform(pars=estpars, transform=transform)
    if (any(is.na(estpars)))
        opt <- list(params=estpars,
                    lik=NA)
    else if (eval.only==TRUE) {
        x <- str_obj(estpars=estpars,
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
                         fn=str_obj,
                         data=obsdata,
                         fixpars=fixpars,
                         parorder=parorder,
                         transform=transform)
        else
            x <- optim(par=estpars,
                       fn=str_obj,
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
str_obj <- function(estpars, data, fixpars, parorder, transform) {
    ## We will give the model the true initial conditions
    y0 <- c(E=100, C=4000, G=0, T=0)
    ## Put the parameters back on the natural scale
    estpars <- par_untransform(estpars, transform)
    ## t_rep lies between 0 and 1 - multiply by max timestep
    estpars["t_rep"] <- estpars["t_rep"] * max(data$time)
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
            dllname="Structured_parasite_model",
            initfunc="initmod",
            events=list(func="event", time=pars["t_rep"])
            )) -> out
    ## if this parameter set produces a simulation error, if it fails
    ## to completely run, or if any of the variables take negative
    ## values, return -Inf for the likelihood.
    if (inherits(out, "try-error"))
        lik <- NA #-Inf
    else if (max(out[,"time"]) < max(data$time))
        lik <- NA #-Inf
    else if (any(is.nan(out)))
        lik <- NA #-Inf
    else if (any(out < 0))
        lik <- NA #-Inf
    else {
        out <- as.data.frame(out)
        ## calculate the log-likelihood of observing these data
        ## assuming only observation error
        sapply(data$time, function(t) out$T[out$time==t]) -> simspores
        dnorm(data$spores, simspores, unname(pars["obs_sd"]), log=T) %>% sum(., na.rm=TRUE) -> lik
    }
    ## return the negative log-likelihood (because optimization methods minimize)
    return(-lik)
}

## trajectory matching function for the unstructured model
traj_match2 <- function(estpars, fixpars, parorder, transform, obsdata, eval.only=FALSE, method="subplex") {
    estpars <- par_transform(pars=estpars, transform=transform)
    if (any(is.na(estpars)))
        opt <- list(params=estpars,
                    lik=NA)
    else if (eval.only==TRUE) {
        x <- unstr_obj(estpars=estpars,
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
                         fn=unstr_obj,
                         data=obsdata,
                         fixpars=fixpars,
                         parorder=parorder,
                         transform=transform)
        else
            x <- optim(par=estpars,
                       fn=unstr_obj,
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
## Objective function to minimize for the unstructured model
unstr_obj <- function(estpars, data, fixpars, parorder, transform) {
    ## Put the parameters back on the natural scale
    estpars <- par_untransform(estpars, transform)
    ## Parasite initial conditions will be estimated
    y0 <- c(E=100, P=unname(estpars["P0"]))
    ## combine the parameters to be estimated and the fixed parameters
    ## into a single vector, with an order specified by parorder
    pars <- c(estpars, fixpars)
    pars[match(parorder, names(pars))] -> pars
    if(any(names(pars)!=parorder)) stop("parameter are in the wrong order")
    holdpars <<- pars
    ## simulate the model at these parameters
    try(ode(y0,
            times=seq(0, max(data$time)),
            func="derivs",
            parms=pars,
            dllname="Unstructured_parasite_model",
            initfunc="initmod"
            )) -> out
    ## if this parameter set produces a simulation error, skip it
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
        ## calculate the log-likelihood of observing these data
        ## assuming only observation error
        sapply(data$time, function(t) out$P[out$time==t]) -> simspores
        dnorm(data$spores, simspores, unname(pars["obs_sd"]), log=T) %>% sum(., na.rm=TRUE) -> lik
    }
    ## return the negative log-likelihood (because optimization methods minimize)
    return(-lik)
}
library(deSolve)
library(magrittr)
library(subplex)
library(pomp)
library(parallel)
library(dplyr)

## Transform the parameters.
## The optimization routines work best if the parameter values are
## unconstrained. Most of the parameters vary between 0 and Inf; these
## can be made unconstrained by simple log-transforming. The parameter
## t_rep, on the other hand, can only vary between 0 and the max time
## in the fitted dataset. Thus, we let t_rep take only values between
## 0 and 1, and then multiply the value of t_rep by the max time to
## get the value of t_rep for the simulations. To put a parameter that
## varies between 0 and 1 on an unconstrained scale, we use the logit
## transform.
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
## Trajectory matching to estimate the parameters of the dynamical model.
## estpars is a named numeric vector giving initial guesses for the
##      parameters whose values will be inferred using trajecory matching
## fixpars is a named numeric vector giving the values of parameters
##      whose values are not estimated
## parorder is a character vector giving the order of parameters
##     expected by the C function that simulates the dynamical model
## transform is a character vector giving the transforms for each of the
##     estimated parameters
## obsdata is a dataframe containing the observed data (times and spore counts)
## eval.only is a logical variable that determines whether the
##     likelihood is to maximized or simply calculated for a set of
##     parameters
## method is the name of the optimizer to use
traj_match <- function(estpars, fixpars, parorder, transform, obsdata, eval.only=FALSE, method="subplex") {
    estpars <- par_transform(pars=estpars, transform=transform)
    if (any(is.na(estpars)))
        opt <- list(params=estpars,
                    lik=NA)
    else if (eval.only==TRUE) {
        x <- str_obj(estpars=estpars,
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
                         fn=str_obj,
                         data=obsdata,
                         fixpars=fixpars,
                         parorder=parorder,
                         transform=transform)
        else
            x <- optim(par=estpars,
                       fn=str_obj,
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
str_obj <- function(estpars, data, fixpars, parorder, transform) {
    ## We will give the model the true initial conditions
    y0 <- c(E=100, C=4000, G=0, T=0)
    ## Put the parameters back on the natural scale
    estpars <- par_untransform(estpars, transform)
    ## t_rep lies between 0 and 1 - multiply by max timestep
    estpars["t_rep"] <- estpars["t_rep"] * max(data$time)
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
            dllname="Structured_parasite_model",
            initfunc="initmod",
            events=list(func="event", time=pars["t_rep"])
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
        ## calculate the log-likelihood of observing these data
        ## assuming only observation error
        sapply(data$time, function(t) out$T[out$time==t]) -> simspores
        dnorm(data$spores, simspores, unname(pars["obs_sd"]), log=T) %>% sum(., na.rm=TRUE) -> lik
    }
    ## return the negative log-likelihood (because optimization methods minimize)
    return(-lik)
}

## trajectory matching function for the unstructured model
traj_match2 <- function(estpars, fixpars, parorder, transform, obsdata, eval.only=FALSE, method="subplex") {
    estpars <- par_transform(pars=estpars, transform=transform)
    if (any(is.na(estpars)))
        opt <- list(params=estpars,
                    lik=NA)
    else if (eval.only==TRUE) {
        x <- unstr_obj(estpars=estpars,
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
                         fn=unstr_obj,
                         data=obsdata,
                         fixpars=fixpars,
                         parorder=parorder,
                         transform=transform)
        else
            x <- optim(par=estpars,
                       fn=unstr_obj,
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
## Objective function to minimize for the unstructured model
unstr_obj <- function(estpars, data, fixpars, parorder, transform) {
    ## Put the parameters back on the natural scale
    estpars <- par_untransform(estpars, transform)
    ## Parasite initial conditions will be estimated
    y0 <- c(E=100, P=unname(estpars["P0"]))
    ## combine the parameters to be estimated and the fixed parameters
    ## into a single vector, with an order specified by parorder
    pars <- c(estpars, fixpars)
    pars[match(parorder, names(pars))] -> pars
    if(any(names(pars)!=parorder)) stop("parameter are in the wrong order")
    holdpars <<- pars
    ## simulate the model at these parameters
    try(ode(y0,
            times=seq(0, max(data$time)),
            func="derivs",
            parms=pars,
            dllname="Unstructured_parasite_model",
            initfunc="initmod"
            )) -> out
    ## if this parameter set produces a simulation error, skip it
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
        ## calculate the log-likelihood of observing these data
        ## assuming only observation error
        sapply(data$time, function(t) out$P[out$time==t]) -> simspores
        dnorm(data$spores, simspores, unname(pars["obs_sd"]), log=T) %>% sum(., na.rm=TRUE) -> lik
    }
    ## return the negative log-likelihood (because optimization methods minimize)
    return(-lik)
}
