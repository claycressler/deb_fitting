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

ndatasets <- 100
dyn.load("Structured_parasite_model.so")
## Simulate data
y0 <- c(E=100, C=4000, G=0, T=0) ## Initial conditions
## Baseline set of parameters - sample around these
pars <- c(theta=10, r=0.1, aC=0.0005, hC=20, aG=0.0001, hG=15, b=20000, m=10000, t_rep=15, obs_sd=100000)
datasets <- vector(mode='list', length=ndatasets)
i <- 1
set.seed(12340987)
while (is.null(datasets[[ndatasets]])) {
    this_pars <- rnorm(length(pars), mean=pars, sd=pars)
    while (any(this_pars < 0))
        this_pars <- rnorm(length(pars), mean=pars, sd=pars)
    names(this_pars) <- names(pars)
    this_pars["theta"] <- 10
    this_pars["t_rep"] <- round(this_pars["t_rep"])
    out <- ode(y0, times=seq(0, 50), func="derivs", parms=this_pars, dllname="Structured_parasite_model", initfunc="initmod", events=list(func="event", time=round(this_pars["t_rep"]))) ## simulated data
    if (all(out >= 0)) {
        data.frame(time=rep(seq(5,50,5), each=10),
                   spores=sapply(out[seq(5,50,5),"T"], function(m) rnorm(10, mean=m, sd=this_pars["obs_sd"])) %>% as.numeric) -> data
        data$spores[data$spores < 0] <- 0
        datasets[[i]] <- list(params=this_pars, data=data)
        i <- i+1
    }
}

str_parorder <- c("theta","r","aC","hC","aG","hG","b","m","t_rep","obs_sd")
unstr_parorder <- c("theta", "r", "aP", "hP", "b", "P0", "obs_sd")
unstr_fixpars <- str_fixpars <- c(theta=10)
str_transform <- rep("log",length(str_parorder)-1)
str_transform[which(str_parorder=="t_rep")-1] <- "logit"
unstr_transform <- rep("log",length(unstr_parorder)-1)

str_parorder <- c("theta","r","aC","hC","aG","hG","b","m","t_rep","obs_sd")
unstr_parorder <- c("theta", "r", "aP", "hP", "b", "P0", "obs_sd")
unstr_fixpars <- str_fixpars <- c(theta=10)
str_transform <- rep("log",length(str_parorder)-1)
str_transform[which(str_parorder=="t_rep")-1] <- "logit"
unstr_transform <- rep("log",length(unstr_parorder)-1)

est_params <- vector(mode='list', length=ndatasets)
for (i in 11:ndatasets) {
    t1 <- Sys.time()
    print(i)
    data <- datasets[[i]]$data

    ## STRUCTURED MODEL
    ## Generate a large number of different initial parameter guesses
    str_box <- cbind(lower=pars[-which(names(pars)=="theta")]/1000, upper=pars[-which(names(pars)=="theta")]*1000)
    str_box[rownames(str_box)=="t_rep",] <- c(0, 1)
    sobolDesign(lower=str_box[,"lower"], upper=str_box[,"upper"], nseq=500000) %>%
        apply(., 1, as.list) %>%
            lapply(., function(l) unlist(l)) -> guesses
    dyn.load("Structured_parasite_model.so")
    mclapply(guesses,
             traj_match,
             fixpars=str_fixpars,
             parorder=str_parorder,
             transform=str_transform,
             obsdata=data,
             eval.only=TRUE,
             method="subplex",
             mc.cores=10) %>%
        lapply(., function(x) x$lik) %>%
            unlist-> guess_lik
    ## Pick the best 1000 of these and use them as starting points for further refinement
    guesses[order(guess_lik)[1:1000]] -> refine
    mclapply(refine,
             traj_match,
             fixpars=str_fixpars,
             parorder=str_parorder,
             transform=str_transform,
             obsdata=data,
             eval.only=FALSE,
             method="Nelder-Mead",
             mc.cores=10) -> refine_lik
    refine_lik %>%
        lapply(., function(l) c(l$params, l$lik)) %>%
            unlist %>%
                matrix(., ncol=10, byrow=TRUE) %>%
                    as.data.frame -> refine_pars
    colnames(refine_pars) = c("r","aC","hC","aG","hG","b","m","t_rep","obs_sd", "lik")
    str_refine_pars <- arrange(refine_pars, lik)

    ## STRUCTURED MODEL
    ## Generate a large number of different initial parameter guesses
    unstr_box <- cbind(lower=c(r=1, aP=1e-4, hP=1, b=100000, P0=1000, obs_sd=1e5)/1000,
                       upper=c(r=1, aP=1e-4, hP=1, b=100000, P0=1000, obs_sd=1e5)*1000)
    sobolDesign(lower=unstr_box[,"lower"], upper=unstr_box[,"upper"], nseq=5000) %>%
        apply(., 1, as.list) %>%
            lapply(., function(l) unlist(l)) -> guesses
    ## compute the likelihood of each of these guesses
    dyn.load("Unstructured_parasite_model.so")
    mclapply(guesses,
             traj_match2,
             fixpars=unstr_fixpars,
             parorder=unstr_parorder,
             transform=unstr_transform,
             obsdata=data,
             eval.only=TRUE,
             method="Nelder-Mead",
             mc.cores=10) %>%
        lapply(., function(l) l$lik) %>%
            unlist -> guess_lik
    ## Pick the best 1000 for further refinement
    guesses[order(guess_lik)[1:1000]] -> refine
    mclapply(refine,
             traj_match2,
             fixpars=unstr_fixpars,
             parorder=unstr_parorder,
             transform=unstr_transform,
             obsdata=data,
             eval.only=FALSE,
             method="Nelder-Mead",
             mc.cores=10) -> refine_lik
    lapply(refine_lik, function(l) c(l$params, l$lik)) %>%
        unlist %>%
            matrix(., ncol=7, byrow=TRUE) %>%
                as.data.frame -> refine_pars
    colnames(refine_pars) = c( "r", "aP", "hP", "b", "P0", "obs_sd", "lik")
    arrange(refine_pars, lik) -> unstr_refine_pars

    est_params[[i]] <- list(str=str_refine_pars, unstr=unstr_refine_pars)
    saveRDS(est_params, file="Comparing_structured_unstructured_model_fits.RDS")
    t2 <- Sys.time()
    print(t2-t1)
}
