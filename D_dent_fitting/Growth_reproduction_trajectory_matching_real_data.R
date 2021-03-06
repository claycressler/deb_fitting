library(deSolve)
library(magrittr)
library(subplex)
library(pomp)
library(parallel)
library(plyr)
library(tidyr)
library(ggplot2)

## I need to think carefully about how to deal with maturation. The DEB model will allow individuals to reproduce immediately. I could treat size at maturity as a parameter to be estimated. Or I could simply say, for each individual, that size at maturity is whatever size reproduction was first observed. I will have to play around with both of these models. I will start with a model where age at maturity has to be estimated along with the other parameters. Thus, size at maturity must be included in the source code for simulating the model.

## rebuild source for this computer
if (is.loaded("tm_deb.so")) dyn.unload("tm_deb.so")
system("rm tm_deb.so")
system("rm tm_deb.o")
system("R CMD SHLIB tm_deb.c")
dyn.load("tm_deb.so")

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

optimizer <- function(estpars, fixpars, parorder, transform, obsdata, Np=100, eval.only=FALSE, type="trajectory_matching", method="subplex") {
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
        opt <- list(params=par_untransform(estpars,transform),
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
                             Np=Np))
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

## Trajectory matching
tm_obj <- function(estpars, data, fixpars, parorder, transform) {
    ## Put the estimated parameters back on the natural scale
    estpars <- par_untransform(estpars, transform)
    ## compute Imax and g based on the estimate of fh
    fixpars["Imax"] <- calc_Imax(unname(estpars["Fh"]))
    fixpars["g"] <- calc_g(unname(estpars["Fh"]))
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

    ## Initial condition: assuming that E/W = rho/v and E+W=ER, then W
    ## + W*rho/v = ER, W(1+rho/v)=ER, W=ER/(1+rho/v)
    y0 <- c(F=unname(pars["F0"]),
            E=0,
            W=unname(pars["ER"]/(1+pars["rho"]/pars["v"])),
            R=0)
    y0["E"] <- unname(y0["W"]*pars["rho"]/pars["v"])

    ## Simulate the system
    try(ode(y0,
            times=seq(0,35,0.1),
            func="derivs",
            parms=pars,
            dllname="tm_deb",
            initfunc="initmod",
            events=list(data=eventdat))) -> out
    if (inherits(out, "try-error"))
        lik <- -Inf
    else if (max(out[,"time"]) < max(data$age))
        lik <- -Inf
    else if (any(is.nan(out)))
        lik <- -Inf
    else if (any(out < 0))
        lik <- -Inf
    ## If no errors, compute the likelihood
    else {
        ## subtract off any reproduction that has happened before the size at maturity was reached
        out[,'R'] <- out[,'R'] - out[max(which(out[,'W'] < pars['Wmat'])),'R']
        out[out[,'R'] < 0,'R'] = 0

        ## extract only the data points that can be compared against the true data
        as.data.frame(out[out[,'time']%in%data$age,]) -> pred

        ## compute the observed weight as Wobs = W + E and compute the observed length prediction as Wobs=xi*Lobs^q; (Wobs/xi)^(1/q)=Lobs
        xi <- 1.8e-3; q <- 3;
        mutate(pred, Wobs=W+E, Lobs=(Wobs/xi)^(1/q)) -> pred

        ## compute the probability of observing the data, given the prediction
        sapply(unique(data$age),
               function(d)
                   c(dnorm(x=data$length[data$age==d],
                           mean=pred$Lobs[pred$time==d],
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

