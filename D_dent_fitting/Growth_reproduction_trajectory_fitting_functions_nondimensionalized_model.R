library(deSolve)
library(magrittr)
library(subplex)
library(pomp)
library(parallel)
library(plyr)
library(tidyr)
library(ggplot2)

if (is.loaded("nondim_deb_Cat.so")) dyn.unload("nondim_deb_Cat.so")
system("rm nondim_deb_Cat.so")
system("R CMD SHLIB nondim_deb_Cat.c")
dyn.load("nondim_deb_Cat.so")

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

calc_Imax <- function(fh) 10616.500976 + 1.129421*fh
calc_g <- function(fh) 1.38413 + 9.839020e-6*fh - 2.738144e-10*fh^2 + 2.648817e-15*fh^3

traj_match <- function(pars, obsdata, events, eval.only=FALSE, method="subplex") {
    ## Some of the parameters have fixed values (v, feeding parameters, length-dry weight parameters)
    fh <- 1e-4/44.5e-9 ## about 2250 cells/ml
    fixpars <- c(Imax=calc_Imax(fh), # fixed based on feeding data fitting and fh value
                 fh=fh, # estimated from fitting growth/reproduction data
                 g=calc_g(fh), # fixed based on feeding data fitting and fh value
                 eps=44.5e-9, # fixed based on measured carbon content of algae
                 V=30, # fixed based on experimental conditions
                 xi=1.8e-3, # fixed based on Hall et al. 2009 length-weight regression
                 q=3, # fixed based on Hall et al. 2009 length-weight regression
                 v=100, # fixed based on the fact that it doesn't affect the fitting
                 F0=1000000/30, # fixed based on experimental conditions
                 L0=0.5, # fixed length at birth
                 Wmat=0.002) # fixed mass at maturity

    ## Some of the dimensionless parameters are estimated
    matrix(pars, nrow=1, byrow=TRUE, dimnames=list(NULL, names(pars))) %>%
        as.data.frame %>%
            with(.,
                 c(w_scalar=EG/K, ## reserve density scalar
                   l_scalar=unname(fixpars["v"])/km, ## length scalar
                   rho=rho, ## dimensionless growth efficiency
                   beta=EG/ER*(1-K)/K*unname(fixpars["xi"])*(unname(fixpars["v"])/km)^unname(fixpars["q"]), ## dimensionless birth rate
                   W0=W0, ## initial reserve density
                   Lobs=Lobs ## observation error for length
                   )
                 ) -> estpars

    ## Transform the parameters to the unconstrained scale
    estpars <- par_transform(pars=estpars, transform=c('log','log','logit','log','log','log'))

    if (any(is.na(estpars)))
        opt <- list(params=pars,
                    lik=NA)
    else if (eval.only==TRUE) {
        x <- obj(estpars=estpars,
                 data=obsdata,
                 fixpars=fixpars,
                 events=events)
        ## Recover the dimensional parameters from the dimensionless ones
        opt <- list(params=pars,
                    lik=x)
    }
    else {
        if (method=="subplex")
            x <- subplex(par=estpars,
                         fn=obj,
                         data=obsdata,
                         fixpars=fixpars,
                         events=events)
        else
            x <- optim(par=estpars,
                       fn=obj,
                       method=method,
                       data=obsdata,
                       fixpars=fixpars,
                       events=events,
                       control=list(maxit=5000))
        opt <- list(params=par_untransform(x$par,c('log','log','logit','log','log','log')),
                    lik=x$value,
                    conv=x$convergence)
    }
    if (eval.only==FALSE) print(opt$lik)
    return(opt)
}
## Objective function to minimize for the structured model
obj <- function(estpars, data, fixpars, events) {
    ## Put the estimated parameters back on the natural scale
    estpars <- par_untransform(estpars, c('log','log','logit','log','log','log'))

    ## State variable scalar multipliers
    scalars <- c(t_scalar=unname(fixpars["fh"]/(fixpars["Imax"]*estpars["l_scalar"]^fixpars["g"])),
                 f_scalar=unname(fixpars["fh"]),
                 estpars["w_scalar"],
                 estpars["l_scalar"])

    ## Dimensionless parameters required for the ODE
    ode_pars <- c(fixpars["g"],
                  fixpars["q"],
                  alpha=unname((fixpars["v"]/estpars["l_scalar"])*scalars["t_scalar"]),
                  estpars["rho"],
                  phi=unname(fixpars["fh"]*fixpars["eps"]*fixpars["V"]/(fixpars["xi"]*estpars["l_scalar"]^fixpars["q"])),
                  estpars["beta"],
                  estpars["w_scalar"])

    ## Dimensionless initial conditions
    y0 <- c(F=unname(fixpars["F0"]/scalars["f_scalar"]),
            W=unname(estpars["W0"]/scalars["w_scalar"]),
            L=unname(fixpars["L0"]/scalars["l_scalar"]),
            R=0)

    ## Dimensionless events
    nondim_events <- events
    nondim_events$time <- round(nondim_events$time/unname(scalars["t_scalar"]))
    nondim_events$value <- nondim_events$value/unname(scalars["f_scalar"])

    ## Dimensionless time steps
    times <- c(0, nondim_events$time)

    ## Simulate the dimensionless model at these parameters
    try(ode(y0,
            times=times,
            func="derivs",
            parms=ode_pars,
            dllname="nondim_deb_Cat",
            initfunc="initmod",
            events=list(data=nondim_events)
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
    else if (any(is.na(out)))
        lik <- -Inf
    else if (any(out < 0))
        lik <- -Inf
    else {
        out <- as.data.frame(out)
        mutate(out,
               time=time*scalars["t_scalar"],
               L=L*scalars["l_scalar"],
               R=R-R[which((fixpars["xi"]*L^fixpars["q"]) < fixpars["Wmat"]) %>% max]) -> out
        out$R[out$R < 0] <- 0
        ## calculate the log-likelihood of observing these data
        ## assuming only observation error
        l_lik <- dnorm(data$length,
                       sapply(data$times,
                              function(t) out$L[floor(out$time)==t]) %>% unlist,
                       unname(estpars["Lobs"]),
                       log=TRUE)
        r_lik <- dpois(data$eggs,
                       sapply(data$times,
                              function(t) out$R[floor(out$time)==t]) %>% unlist,
                       log=TRUE)
        if (length(l_lik) < 1 | length(r_lik) < 1) ## likelihood cannot be evaluated
            lik <- -Inf
        else
            sum(l_lik+r_lik, na.rm=TRUE) -> lik
    }
    ## return the negative log-likelihood (because optimization methods minimize)
    return(-lik)
}


