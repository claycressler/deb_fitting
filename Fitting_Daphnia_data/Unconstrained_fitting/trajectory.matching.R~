require(deSolve)


## I have written my own code for trajectory matching, as it is not
## possible inside of pomp() for batch transfer conditions. This is
## because the food state variable must be updated discontinuously,
## which is disallowed in lsoda, which is the ode integrator used by
## pomp. I have written up the batch transfer equations in C-code that
## can be compiled and linked. This file includes the necessary code
## for discontinuous changes in the food variable.
## Here is my function. Its arguments are very similar to those of
## traj.match in the pomp package.
## "data" is a dataframe of observed data to fit the model to.
## "startpars" gives starting parameter estimates on the natural
##   scale; this includes the initial conditions, which we want to be
##   able to estimate as well
## "estpars" gives the names of the parameters to estimate.
## "method" is the optimization method - choices are any of the
##   methods used by 'optim'
trajectory.matching <- function(data, pars, transform=NA,
                                est, method='subplex',
                                odefunc, initfunc, dllname, 
                                eval.only=FALSE, maxit=500, nondim=TRUE, ...) {
    if (!is.data.frame(data)) stop('`data` must be of class `dataframe`')
    if (!is.numeric(pars)) {
        p <- as.numeric(pars)
        names(p) <- names(pars)
        p -> pars
    }

    ## Load the dynamic library containing the function that numerically
    ## integrates the batch transfer ODE modeldf
    dyn.load(paste(dllname,'so',sep='.'))

    ## Record the order of the parameter values in startpars - these
    ## must match the order expected by the C-function
    name.order <- names(pars)

    if (length(transform) > 1)
        ## Transform parameters to the estimation scale using the
        ## transformation specified by "transform"
        pars <- est.scale.params(pars, transform)

    ## Extract the initial values of the parameters to be optimized
    ## over
    est.pars <- pars[sapply(1:length(est),function(i) match(est[i],names(pars)))]

    ## Remove these parameters from the full parameter vector to avoid
    ## any potential conflicts (overwriting of parameter values) in
    ## the objective function
    fixed.pars <- pars[-(sapply(1:length(est),function(i) match(est[i],names(pars))))]

    if (eval.only==TRUE) {
        x <- obj.func(estpars=est.pars, data=data, fixedpars=fixed.pars, parorder=name.order, transform=transform, odefn=odefunc, initfn=initfunc, dll=dllname, nondim=nondim, ...)
        opt <- list(par=nat.scale.params(pars,transform), lik=x)
        return(opt)
    }
    else {
        if (method=='subplex')
            x <- subplex(par=est.pars, fn=obj.func, data=data, fixedpars=fixed.pars, parorder=name.order, transform=transform, odefn=odefunc, initfn=initfunc, dll=dllname, nondim=nondim, control=list(maxit=maxit), ...)
        else
            x <- optim(est.pars, fn=obj.func, method=method, data=data, fixedpars=fixed.pars, parorder=name.order, transform=transform, odefn=odefunc, initfn=initfunc, dll=dllname, nondim=nondim, control=list(maxit=maxit), ...)

        parameters <- sort.params(fixed.pars, x$par, name.order)
        opt <- list(par=nat.scale.params(parameters, transform),
                    lik=x$value,
                    conv=x$convergence)

        return(opt)
    }
}

## Functions for rescaling the parameters according to "transform"
est.scale.params <- function(p, transform) {
    if (length(p)!=length(transform))
        stop('User must specify transformation rules for all parameters')
    for (i in 1:length(p)) {
        if (transform[i]=='log')
            p[i] <- log(p[i])
        else if (transform[i]=='logit')
            p[i] <- log(p[i]/(1-p[i]))
    }
    return(p)
}

nat.scale.params <- function(p, transform) {
    if (length(p)!=length(transform))
        stop('User must specify transformation rules for all parameters')
    for (i in 1:length(p)) {
        if (transform[i]=='log')
            p[i] <- exp(p[i])
        else if (transform[i]=='logit')
            p[i] <- exp(p[i])/(1+exp(p[i]))
    }
    return(p)
}

sort.params <- function(parset1, parset2, order) {
    sorted.params <- vector()
    for (i in 1:length(order)) {
        if (!is.na(match(order[i],names(parset1))))
            sorted.params <- c(sorted.params,
                               parset1[match(order[i],names(parset1))])
        else
            sorted.params <- c(sorted.params,
                               parset2[match(order[i],names(parset2))])
    }
    return(sorted.params)
}

obj.func <- function(estpars, data, fixedpars, parorder, transform, odefn, initfn, dll, nondim,...) {

    ## Merge estpars and fixedpars, matching parorder
    params <- sort.params(fixedpars, estpars, parorder)

    ## Convert the parameters to the natural scale; kappa must be the
    ## first parameter listed
    if (length(transform) > 1)
        params <- nat.scale.params(params, transform)

    ## Extract the initial value parameters
    y0 <- params[grep('.0',names(params))]
    names(y0) <- gsub('.0','',names(y0))

    ## The time points to return the deterministic trajectory match
    ## the time points of observation in the dataset to be fitted If
    ## the model is nondimensional, then the length of time depends on
    ## the value of the time scalar.
    times <- data$age
    if (nondim==TRUE) times <- times*unname(params['time.scalar'])

    ## Calculate the deterministic trajectory
    out <- try(ode(y0, c(0,times), func=odefn, parms=params, dllname=dll, initfunc=initfn,...))
    
    if (inherits(out, 'try-error')) {
        warning('Error in ODE integrator')
        lik <- NA
      }
    else {
        ## Remove the initial condition row from out
        out <- as.data.frame(out)
        out <- out[-1,]

        ## Redimensionalize the output, if necessary
        if (length(transform) > 1) {
            L <- out$Lmeas*unname(params['L.scalar'])
            R <- out$Rm-unname(params['beta'])
            R[which(R <= 0)] <- 0
            R <- R*unname(params['Rm.scalar'])
          
        }

        ## Calculate the likelihood of the output given the input
        L.sd <- unname(params['L.sd'])
        R.sd <- unname(params['R.sd'])
        if (nrow(out) != nrow(data))
            lik <- NA
        else
            lik <- sum(dnorm(data$Lobs, L, L.sd, log=T)
                       + dnorm(data$Robs, R, R.sd, log=T), na.rm=T)
    }
    return(-lik)
}
