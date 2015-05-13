require(deSolve)

## Load the dynamic library containing the function that numerically
## integrates the batch transfer ODE modeldf
dyn.load("deb_batch.so")
dyn.load("deb_constant.so")

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
## "fixed.pars" is a named vector of parameter values for all fixed
##     (non-estimated) parameters on the natural scale.
## "estd.pars" is a named vecctor of initial parameter estimates on
##      the natural scale
##  (Note that both "fixed.pars" and "estd.pars" may contain initial
##  conditions.)
## "method" is the optimization method - choices are any of the
##      methods used by 'optim'
## "eval.only" determines whether trajectory matching occurs or
##      whether the likelihood is evaluated at the starting parameters
##      and returned
## "maxit" is the number of optimizer iterations to run
## "model" determines whether the fit model is for batch culture or
##      constant ingestion conditions

trajectory.matching <- function(data, fixed.pars, estd.pars, method='subplex', eval.only=FALSE, maxit=500, model) {
    if (!is.data.frame(data)) stop('`data` must be of class `dataframe`')

    ## The order of the parameters must match the order expected by
    ## the C-function
    name.order <- c('K','km','eG','eR','v','Rmbar','Imax','Fh','eA','E.0','L.0','R.0','F.0','L.sd')

    ## Transform parameters to the estimation scale - log-transform
    ## all but kappa and eA, which are logit-transformed
    estd.pars <- est.scale.params(estd.pars)
    fixed.pars <- est.scale.params(fixed.pars)

    if (eval.only==TRUE) {
        x <- obj.func(estpars=estd.pars, data=data, fixedpars=fixed.pars, parorder=name.order, model=model)
        opt <- list(par=nat.scale.params(sort.params(fixed.pars, estd.pars, name.order)), lik=x)
        return(opt)
    }
    else {
        if (method=='subplex')
            x <- subplex(par=estd.pars, fn=obj.func, data=data, fixedpars=fixed.pars, parorder=name.order, model=model, control=list(maxit=10000))
        else
            x <- try(optim(estd.pars, fn=obj.func, method='Nelder-Mead', data=data, fixedpars=fixed.pars, parorder=name.order, model=model, control=list(maxit=maxit)))
        if (inherits(x, 'try-error'))
          opt <- list(par=nat.scale.params(sort.params(fixed.pars, estd.pars, name.order)),
                      lik=NA,
                      conv=NA)
        else
          opt <- list(par=nat.scale.params(sort.params(fixed.pars, x$par, name.order)),
                      lik=x$value,
                      conv=x$convergence)
        return(opt)
    }
}

est.scale.params <- function(p) {
    p <- log(p)
    ## logit transform K and eA
    if ('K'%in%names(p))
        p['K'] <- p['K']-log(1-exp(p['K']))
    if ('eA'%in%names(p))
        p['eA'] <- p['eA']-log(1-exp(p['eA']))
    return(p)
}

nat.scale.params <- function(p) {
    p <- exp(p)
    ## expit transform K and eA
    if ('K'%in%names(p))
        p['K'] <- p['K']/(1+p['K'])
    if ('eA'%in%names(p))
        p['eA'] <- p['eA']/(1+p['eA'])
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

obj.func <- function(estpars, data, fixedpars, parorder, model) {
    ## Merge estpars and fixedpars, matching parorder
    params <- sort.params(fixedpars, estpars, parorder)

    ## Convert the parameters to the natural scale; kappa must be the
    ## first parameter listed
    params <- nat.scale.params(params)

    ## The time points to return the deterministic trajectory match
    ## the time points of observation in the dataset to be fitted
    times <- data$age

    ## Extract the initial value parameters
    y0 <- params[grep('.0',names(params))]
    names(y0) <- gsub('.0','',names(y0))

    ## Calculate the deterministic trajectory
    if (model=='batch')
        out <- try(ode(y0, c(0,times), func="derivs", parms=params, dllname="deb_batch", initfunc="initmod", events=list(func="event",time=times),maxsteps=5000))
    if (model=='constant') {
        ## F0 is modified to F0/t_TI, where t_TI is the length of the
        ## transfer interval
        params['F.0'] <- params['F.0']/diff(times)[1]
        y0['F'] <- y0['F']/diff(times)[1]
        out <- try(ode(y0, c(0,times), func="derivs", parms=params, dllname="deb_constant", initfunc="initmod", maxsteps=5000))
    }
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
