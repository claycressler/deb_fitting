library(LambertW)
library(magrittr)
library(subplex)
library(pomp)
library(parallel)
library(plyr)

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
## Calculate the algal concentration at any time point, for given value of the half-saturation constant
calc_Ft <- function(a, Z, h, Imax, g, F0, L, fh, t)
    fh * W(exp(F0/fh - exp(-a*Z/L^h)*Imax*L^g*t/fh) * F0 / fh)


traj_match <- function(estpars, fixpars, transform, obsdata, eval.only=FALSE, method="subplex") {
    estpars <- par_transform(pars=estpars, transform=transform)
    if (any(is.na(estpars)))
        opt <- list(params=estpars,
                    lik=NA)
    else if (eval.only==TRUE) {
        x <- obj(estpars=estpars,
                 data=obsdata,
                 fixpars=fixpars,
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
                         transform=transform)
        else
            x <- optim(par=estpars,
                       fn=obj,
                       method=method,
                       data=obsdata,
                       fixpars=fixpars,
                       transform=transform,
                       control=list(maxit=5000))
        opt <- list(params=par_untransform(x$par,transform),
                    lik=x$value,
                    conv=x$convergence)
    }
    return(opt)
}
## Objective function to minimize for the feeding model
obj <- function(estpars, data, fixpars, transform) {
    ## Put the estimated parameters back on the natural scale
    estpars <- par_untransform(estpars, transform)
    pars <- c(estpars, fixpars)
    ## compute the expected amount of food left
    mutate(data,
           E.F0=pars["F0"],
           E.Ft=calc_Ft(pars["a"],Z,pars["h"],pars["Imax"],pars["g"],pars["F0"],length,pars["fh"],t)) -> ndata
    with(ndata,
         dnorm(F0, mean=E.F0, sd=pars["Fobs"], log=TRUE) +
             dnorm(Ft, mean=E.Ft, sd=pars["Fobs"], log=TRUE)
         ) %>% sum -> lik
    ## return the -log(L) (b/c optimization methods minimize)
    return(-lik)
}
