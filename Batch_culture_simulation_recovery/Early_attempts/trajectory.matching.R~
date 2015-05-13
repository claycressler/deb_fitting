require(deSolve)

## Load the dynamic library containing the function that numerically
## integrates the batch transfer ODE modeldf
dyn.load("my_deb_batch_constantK.so")

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
trajectory.matching <- function(data, pars, est, method='subplex', eval.only=FALSE, maxit=500) {
    if (!is.data.frame(data)) stop('`data` must be of class `dataframe`')
    if (!is.numeric(pars)) {
        p <- as.numeric(pars)
        names(p) <- names(pars)
        p -> pars
    }

    ## Record the order of the parameter values in startpars - these
    ## must match the order expected by the C-function
    name.order <- names(pars)

    ## Transform parameters to the estimation scale - log-transform
    ## all but kappa, which is logit-transformed
    pars <- est.scale.params(pars)

    ## Extract the initial values of the parameters to be optimized
    ## over
    est.pars <- pars[sapply(1:length(est),function(i) grep(est[i],names(pars)))]

    ## Remove these parameters from the full parameter vector to avoid
    ## any potential conflicts (overwriting of parameter values) in
    ## the objective function
    fixed.pars <- pars[-sapply(1:length(est),function(i) grep(est[i],names(pars)))]

    if (eval.only==TRUE) {
        x <- obj.func(estpars=est.pars, data=data, fixedpars=fixed.pars, parorder=name.order)
        opt <- list(par=nat.scale.params(pars), lik=x)
        return(opt)
    }
    else {
        if (method=='subplex')
            x <- subplex(par=est.pars, fn=obj.func, data=data, fixedpars=fixed.pars, parorder=name.order, control=list(maxit=10000))
        else
            x <- optim(est.pars, fn=obj.func, method=method, data=data, fixedpars=fixed.pars, parorder=name.order, control=list(maxit=maxit))

        opt <- list(par=nat.scale.params(sort.params(fixed.pars, x$par, name.order)),
                    lik=x$value,
                    conv=x$convergence)

        return(opt)
    }
}

est.scale.params <- function(p) {
    p['kappa'] <- log(p['kappa']/(1-p['kappa']))
    p['eA'] <- log(p['eA']/(1-p['eA']))
    p[-which(names(p)%in%c('kappa','eA'))] <- log(p[-which(names(p)%in%c('kappa','eA'))])
    return(p)
}

nat.scale.params <- function(p) {
    p['kappa'] <- exp(p['kappa'])/(1+exp(p['kappa']))
    p['eA'] <- exp(p['eA'])/(1+exp(p['eA']))
    p[-which(names(p)%in%c('kappa','eA'))] <- exp(p[-which(names(p)%in%c('kappa','eA'))])
    return(p)
}

sort.params <- function(parset1, parset2, order) {
    sorted.params <- vector()
    for (i in 1:length(order)) {
        if (length(grep(order[i],names(parset1)))>0)
            sorted.params <- c(sorted.params,
                               parset1[grep(order[i],names(parset1))])
        else
            sorted.params <- c(sorted.params,
                               parset2[grep(order[i],names(parset2))])
    }
    return(sorted.params)
}

obj.func <- function(estpars, data, fixedpars, parorder) {
    ## Merge estpars and fixedpars, matching parorder
    params <- sort.params(fixedpars, estpars, parorder)

    ## Convert the parameters to the natural scale; kappa must be the
    ## first parameter listed
    params <- nat.scale.params(params)

    ## Extract the initial value parameters
    y0 <- params[grep('.0',names(params))]
    names(y0) <- gsub('.0','',names(y0))

    ## The time points to return the deterministic trajectory match
    ## the time points of observation in the dataset to be fitted
    times <- data$time

    ## Calculate the deterministic trajectory
    out <- try(ode(y0, c(0,times), func="derivs", parms=params, dllname="my_deb_batch_constantK", initfunc="initmod", events=list(func="event",time=times),maxsteps=5000))

    if (inherits(out, 'try-error'))
        stop('Error in ODE integrator')
    else {
        ## Remove the initial condition row from out
        out <- as.data.frame(out)
        out <- out[-1,]
        ## Calculate the likelihood of the output given the input
        L.sd <- unname(params[grep('L.sd', names(params))])
        if (nrow(out) != nrow(data))
            lik <- NA
        else
            lik <- sum(dnorm(data$Lobs, out$L, L.sd, log=T)
                       + dpois(data$Robs, out$R, log=T))
    }
    return(-2*lik)
}
