library(deSolve)
library(magrittr)
library(subplex)
library(pomp)
library(parallel)
library(plyr)
library(tidyr)
library(ggplot2)

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
            W=unname(5.39e-5/(1+pars["rho"]/pars["v"])),
            P1=unname(pars["P0"]),
            P2=0)
    y0["E"] <- unname(y0["W"]*pars["rho"]/pars["v"])

    ## Simulate the system
    try(dede(y0,
            times=seq(0,35,0.1),
            func=delay_model_1,
            parms=pars,
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
        sapply(unique(data$age),
               function(d)
                   c(dnorm(x=data$length[data$age==d],
                           mean=pred$Lobs[pred$time==d],
                           sd=pars["Lobs"],
                           log=TRUE) %>% sum,
                     dnbinom(x=data$spores[data$age==d],
                             mu=pred$P2[pred$time==d],
                             size=pars["Pobs"],
                             log=TRUE) %>% sum
                     ) %>% sum
               ) %>% sum -> lik
    }
    return(-lik)
}

delay_model_1 <- function(t, y, parms) {
    rho <- parms["rho"] ## host assimilation efficiency
    K <- parms["K"] ## host growth allocation
    km <- parms["km"] ## host maintenance rate
    v <- parms["v"] ## host energy conductance
    F0 <- parms["F0"] ## food addition
    Lerr <- parms["Lobs"] ## error in the observation of length
    aP <- parms["aP"] ## parasite's attack rate
    eP <- parms["eP"] ## parasite's conversion efficiency
    kP <- parms["kP"] ## parasites per unit body mass
    tau <- parms["tau"] ## development time for parasite
    P0 <- parms["P0"] ## initial parasite abundance
    Perr <- parms["Pobs"]

    eps <- 0.00000000816536; ## carbon content of algae from Meg's data
    V <- 30; ## volume of the container
    xi <- 0.0018; ## length-weight regression coefficient from Spencer
    q <- 3; ## length-weight regression exponent from Spencer

    ## foraging-dependent parameters
    Fh <- 12000; ## half-saturation constant
    Imax <- 12162; ## ingestion rate
    g <- 1.467; ## size-dependence of ingestion
    a <- 3.2e-5; ## spore coefficient
    h <- 3.57; ## size -dependence of spore dependence

    ## system state
    F <- y[1]
    E <- y[2]
    W <- y[3]
    P1 <- y[4]
    P2 <- y[5]
    L <- W^(1/3)
    Lobs <- (W/xi)^(1/q)

    ## Ingestion
    ing <- Imax * F/(Fh+F) * Lobs^g * exp(-a * P2/Lobs^h)
    ## Mobilization
    pc <- E * (v/L + km) / (1 + K*E/W)

    ## state tau units in the past
    dF <- -ing
    dE <- rho*eps*V*ing - pc - aP*E*P1
    dW <- K*pc - km*W
    if (t < tau) {
        dP1 <- eP*aP*E*P1*(1-P1/(kP*W))
        dP2 <- 0
    }
    else {
        lag <- lagvalue(t-tau)
        dP1 <- eP*aP*E*P1*(1-P1/(kP*W))-eP*aP*lag[2]*lag[4]*(1-lag[4]/(kP*lag[3])) ## maturation
        dP2 <- eP*aP*lag[2]*lag[4]*(1-lag[4]/(kP*lag[3]))
    }

    list(c(dF,dE,dW,dP1,dP2))
}
