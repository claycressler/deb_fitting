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
calc_Ft <- function(pars, L, fh, t) {
    Imax <- pars['Imax']
    g <- pars['g']
    F0 <- pars['F0']

    fh * W(exp(F0/fh - Imax*L^g*t/fh) * F0 / fh)
}

traj_match <- function(estpars, fh, transform, obsdata, eval.only=FALSE, method="subplex") {
    estpars <- par_transform(pars=estpars, transform=transform)
    if (any(is.na(estpars)))
        opt <- list(params=estpars,
                    lik=NA)
    else if (eval.only==TRUE) {
        x <- obj(estpars=estpars,
                 data=obsdata,
                 fh=fh,
                 transform=transform)
        opt <- list(params=estpars,
                    lik=x)
    }
    else {
        if (method=="subplex")
            x <- subplex(par=estpars,
                         fn=obj,
                         data=obsdata,
                         fh=fh,
                         transform=transform)
        else
            x <- optim(par=estpars,
                       fn=obj,
                       method=method,
                       data=obsdata,
                       fh=fh,
                       transform=transform,
                       control=list(maxit=5000))
        opt <- list(params=par_untransform(x$par,transform),
                    lik=x$value,
                    conv=x$convergence)
    }
    return(opt)
}
## Objective function to minimize for the feeding model
obj <- function(estpars, data, fh, transform) {
    ## Put the estimated parameters back on the natural scale
    estpars <- par_untransform(estpars, transform)
    ## compute the expected amount of food left
    mutate(data,
           E.F0=estpars["F0"],
           E.Ft=calc_Ft(estpars,length,fh,t)) -> ndata
    with(ndata,
         dnorm(F0, mean=E.F0, sd=estpars["Fobs"], log=TRUE) +
             dnorm(Ft, mean=E.Ft, sd=estpars["Fobs"], log=TRUE)
         ) %>% sum -> lik
    ## return the -log(L) (b/c optimization methods minimize)
    return(-lik)
}

data <- read.csv("Cat_data/feeding2.csv")
estpars <- c(Imax=1e5, g=2, F0=1e7, Fobs=1e4)
transform <- rep("log",4)
fh <- 2.5e6
## TESTING
obj(par_transform(estpars, transform), data, fh, transform)
traj_match(estpars, fh, transform, data, eval.onl=TRUE)
mutate(data,
       E.F0=estpars["F0"],
       E.Ft=calc_Ft(estpars,length,fh,t)) %>%
    with(.,
         dnorm(F0, mean=E.F0, sd=estpars["Fobs"], log=TRUE) +
             dnorm(Ft, mean=E.Ft, sd=estpars["Fobs"], log=TRUE)
         ) %>% sum

data <- subset(data, infected==0)
## Ranges for the estimated parameters
box <- cbind(lower=c(Imax=1e2, g=1, F0=5000, Fobs=1e2),
             upper=c(Imax=1e5, g=4, F0=15000, Fobs=1e3))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=100000) %>%
    apply(., 1, as.list) %>%
        lapply(., unlist) -> guesses

## Pick a huge range of potential fh values, and perform the fitting
## estimating the other parameters
fh_vals <- seq(100, 10000, 100)
est_params <- vector(mode='list', length=length(fh_vals))
for (i in 1:length(fh_vals)) {
    fh <- fh_vals[i]
    print(fh)
    t1 <- Sys.time()
    mclapply(guesses,
             traj_match,
             fh,
             transform,
             data,
             eval.only=TRUE,
             mc.cores=5) %>%
        lapply(., function(x) x$lik) %>%
            unlist -> guess_lik
    t2 <- Sys.time()
    print(paste("Guess lik calc", t2-t1))
    guesses[order(guess_lik)[1:500]] -> refine
    t1 <- Sys.time()
    mclapply(refine,
             traj_match,
             fh,
             transform,
             data,
             eval.only=FALSE,
             method="Nelder-Mead",
             mc.cores=5) %>%
        lapply(., unlist) %>%
            unlist %>%
                matrix(., ncol=6, byrow=TRUE) %>%
                    as.data.frame -> refine_lik
    colnames(refine_lik) <- c(names(estpars),'lik','conv')
    print(min(refine_lik$lik))
    t2 <- Sys.time()
    print(paste("Refine lik calc", t2-t1))
    est_params[[i]] <- refine_lik
    saveRDS(est_params, "Feeding_model_fitting_2.RDS")
}

