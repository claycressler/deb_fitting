## This differs from Growth_reproduction_trajectory_matching_real_data.R in that it assumes that reproduction only starts once some level of maturity has been reached.
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
if (is.loaded("tm_deb_2.so")) dyn.unload("tm_deb_2.so")
system("rm tm_deb_2.so")
system("rm tm_deb_2.o")
system("R CMD SHLIB tm_deb_2.c")
dyn.load("tm_deb_2.so")

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

optimizer <- function(estpars, fixpars, parorder, transform, obsdata, Np=100, eval.only=FALSE, errmodel=c("normal","normal_cv","negbinom"), type="trajectory_matching", method="subplex") {
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
                        transform=transform,
                        err=errmodel)
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
                                 transform=transform,
                                 err=errmodel))
            else
                x <- try(optim(par=estpars,
                               fn=tm_obj,
                               method=method,
                               data=obsdata,
                               fixpars=fixpars,
                               parorder=parorder,
                               transform=transform,
                               err=errmodel,
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
tm_obj <- function(estpars, data, fixpars, parorder, transform, err) {
     ## Put the estimated parameters back on the natural scale
    estpars <- par_untransform(estpars, transform)
    ## combine the parameters to be estimated and the fixed parameters
    ## into a single vector, with an order specified by parorder
    pars <- c(estpars, fixpars)
    pars[match(parorder, names(pars))] -> pars
    if(any(names(pars)!=parorder)) stop("parameter are in the wrong order")

    ## compute Imax and g based on the estimate of fh
    pars["Imax"] <- calc_Imax(unname(pars["Fh"]))
    pars["g"] <- calc_g(unname(pars["Fh"]))
    print(pars)

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
            M=0,
            R=0)
    y0["E"] <- unname(y0["W"]*pars["rho"]/pars["v"])

    ## Simulate the system
    try(ode(y0,
            times=seq(0,35,0.1),
            func="derivs",
            parms=pars,
            dllname="tm_deb_2",
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
        ## extract only the data points that can be compared against the true data
        as.data.frame(out[out[,'time']%in%data$age,]) -> pred

        ## compute the observed weight as Wobs = W + E and compute the observed length prediction as Wobs=xi*Lobs^q; (Wobs/xi)^(1/q)=Lobs
        xi <- 1.8e-3; q <- 3;
        mutate(pred, Wobs=W+E, Lobs=(Wobs/xi)^(1/q)) -> pred

        ## compute the probability of observing the data, given the prediction
        if (err=="normal")
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
        else if (err=="normal_cv") {
            ## Have to be more careful here: if the model predicts 0 reproduction, that makes the sd of the normal ddistribution 0 as well, which means that the log-likelihood will be Inf. So compute the reproduction likelihood carefully - possibly by changing all zeros to some small non-zero number.
            sapply(unique(data$age),
                   function(d)
                       c(dnorm(x=data$length[data$age==d],
                               mean=pred$Lobs[pred$time==d],
                               sd=pars["Lobs"],
                               log=TRUE) %>% sum,
                         dnorm(x=data$eggs[data$age==d],
                               mean=(pred$R[pred$time==d]+0.001),
                               sd=(pred$R[pred$time==d]+0.001)*pars["Robs"],
                               log=TRUE) %>% sum
                         ) %>% sum
                   ) %>% sum -> lik
        }
        else
            sapply(unique(data$age),
                   function(d)
                       c(dnorm(x=data$length[data$age==d],
                               mean=pred$Lobs[pred$time==d],
                               sd=pars["Lobs"],
                               log=TRUE) %>% sum,
                         dnbinom(x=data$eggs[data$age==d],
                                 mu=pred$R[pred$time==d],
                                 size=pars["Robs"],
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
    ## combine the parameters to be estimated and the fixed parameters
    ## into a single vector, with an order specified by parorder
    pars <- c(estpars, fixpars)
    pars[match(parorder, names(pars))] -> pars
    if(any(names(pars)!=parorder)) stop("parameter are in the wrong order")

    ## compute Imax and g based on the estimate of fh
    pars["Imax"] <- calc_Imax(unname(pars["Fh"]))
    pars["g"] <- calc_g(unname(pars["Fh"]))

    ## Generic set of food addition events
    eventdat <- data.frame(var="F",
                           time=1:35,
                           value=unname(pars['F0']),
                           method=rep(c(rep("add",4),"rep"),max(data$age)/5))

    ## Generate Np sets of initial conditions to initilize both the filtering distribution and the prediction distribution
    ## Initial condition: assuming that E/W = rho/v and E+W=ER, then W
    ## + W*rho/v = ER, W(1+rho/v)=ER, W=ER/(1+rho/v)
    ## error in food addition comes from fitting of feeding data
    Ferr <- 2000
    data.frame(T=0,
               F=rnorm(Np, mean=unname(pars["F0"]), sd=Ferr),
               E=0,
               W=unname(pars["ER"]/(1+pars["rho"]/pars["v"])),
               M=0,
               R=0) %>%
                   mutate(., E = W*pars['rho']/pars['v']) -> x.F -> x.P

    ## timepoints when the likelihood should be evaluated
    times <- c(0, data$age %>% unique)
    tstep <- 1
    lik <- 0
    while (tstep < length(times)) {
        print(tail(x.F))
        ## For each of the Np particles
        for (i in 1:Np) {
            ## generate the food addition events
            events <- subset(eventdat, time >= times[tstep] & time < times[tstep+1])
            events$value <- rnorm(nrow(events),
                                  mean=events$value,
                                  sd=Ferr)
            # obtain a sample of points from the prediction distribution by simulating the model forward
            try(ode(x.F[i,2:6] %>% unlist,
                    times=seq(times[tstep], times[tstep+1]),
                    func="derivs",
                    parms=pars,
                    dllname="tm_deb_2",
                    initfunc="initmod",
                    events=list(data=events))) -> out
            if (inherits(out, "try-error") ||
                max(out[,"time"]) < times[tstep+1] ||
                any(is.nan(out)) ||
                any(out < -1e-4))
                x.P[i,] <- rep(NA,6)
            else tail(out,1) -> x.P[i,]
        }

        ## determine the weights by computing the probability of observing the data, given the points in the prediction distribution
        xi <- 1.8e-3; q <- 3
        sapply((x.P$W/xi)^(1/q),
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
