source("Growth_reproduction_trajectory_fitting_stochastic_functions.R")

pars <- c(Imax=22500, fh=10000, g=1.45, rho=0.1, eps=44.5e-9, V=30, F0=1000000/30, xi=2.62e-3, q=2.4, K=0.3, km=0.15, ER=1.51e-3, v=10, Lobs=0.1, Robs=1, Ferr=5000)
y0 <- c(F=1000000/30, E=0.00025, W=0.00025, R=0)

## days when observations will take place
times <- c(5,10,12,15,18,25,30,35)
## number of reps per observation day
inds <- 1:12

set.seed(101)
out <- vector(mode='list', length=length(times)*length(inds))
for (i in 1:(length(days)*length(inds))) {
    ## feeding schedule - amount of food added each time is stochastic
    eventdat <- data.frame(var="F",
                           time=1:35,
                           value=rnorm(35,
                               mean=unname(pars["F0"]),
                               sd=unname(pars["Ferr"])),
                           method=rep(c(rep("add",4),"rep"),max(times)/5))
    ode(y0, times=0:35, func="derivs", parms=pars, dllname="debStochEnv", initfunc="initmod", events=list(data=eventdat)) %>% as.data.frame -> out[[i]]
}

## Simulate observations
set.seed(101)
data = expand.grid(ind=inds, age=times)
data$length <- rnorm(nrow(data),
                     mean=(sapply(1:nrow(data), function(i) out[[i]][data[i,'age']+1,4])/pars['xi'])^(1/pars['q']),
                     sd=pars['Lobs'])
data$eggs <- rnorm(nrow(data),
                   mean=sapply(1:nrow(data), function(i) out[[i]][data[i,'age']+1,5]),
                   sd=pars['Robs'])

## Implement a particle filter to compute the likelihood of a particular parameter set.
eventdat <- data.frame(var="F",
                       time=1:35,
                       value=1e6/30,
                       method=rep(c(rep("add",4),"rep"),max(times)/5))

times <- c(0,5,10,12,15,18,25,30,35)
## Generate Np sets of initial conditions for the state variables
x <- vector(mode='list', length=Np)
for (i in 1:length(x)) {
    ## There is variance in size at birth that is distinct from the error in the observation of length, which I will ignore for now
    x[[i]] <- array(NA, dim=c(length(times), 5))
    x[[i]][1,] <- c(0,
                    rnorm(1, mean=1e6/30, sd=pars['Ferr']),
                    0.00025,
                    0.00025,
                    0)
    colnames(x[[i]]) <- c('T','F','E','W','R')
}
tstep <- 1

## Simulate to the next timestep
while (tstep < length(times)) {
    ## For each of the Np particles
    for (i in 1:Np) {
        ## generate the food addition events
        events <- subset(eventdat, time >= times[tstep] & time < times[tstep+1])
        events$value <- rnorm(nrow(events),
                              mean=events$value,
                              sd=pars['Ferr'])
        ## iterate the state variables from the current time to the next timestep and save the result
        ode(x[[i]][tstep,2:5],
            times=seq(times[tstep], times[tstep+1], 0.001),
            func="derivs",
            parms=pars,
            dllname="debStochEnv",
            initfunc="initmod",
            events=list(data=events)) %>%
                tail(1) -> x[[i]][tstep+1,]
    }
    ## determine the weights by computing the likelihood of observing the data up to and including the current time






box <- cbind(lower=c(fh=2000, rho=0, K=0, km=0.001, Lobs=0.001, Robs=0.01, Ferr=100),
             upper=c(fh=20000, rho=1, K=1, km=1, Lobs=2, Robs=10, Ferr=10000))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=100000) %>%
    apply(., 1, as.list) %>%
        lapply(., unlist) -> guesses

fixpars <- c(Imax=22500, g=1.45, eps=44.5e-9, V=30, F0=1000000/30, xi=2.62e-3, q=2.4, ER=1.51e-3, v=100)
estpars <- guesses[[1]]
transform <- c("log", rep("logit",2), rep("log",4))
parorder <- c("Imax","fh","g","rho","eps","V","F0","xi","q","K","km","ER","v","Lobs","Robs","Ferr")



