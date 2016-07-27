source("Growth_reproduction_trajectory_fitting_stochastic_functions.R")

pars <- c(Imax=22500, fh=10000, g=1.45, rho=0.1, eps=44.5e-9, V=30, F0=1000000/30, xi=2.62e-3, q=2.4, K=0.3, km=0.15, ER=1.51e-3, v=10, Lobs=0.1, Robs=1, Ferr=5000)
y0 <- c(F=1000000/30, E=0.00025, W=0.00025, R=0)

## days when observations will take place
times <- c(5,10,12,15,18,25,30,35)
## number of reps per observation day
inds <- 1:12

set.seed(101)
out <- vector(mode='list', length=length(times)*length(inds))
for (i in 1:(length(times)*length(inds))) {
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
## Generate Np sets of initial conditions from the filtering distribution
Np <- 100
data.frame(T=0,
           F=rnorm(Np, mean=1e6/30, sd=pars['Ferr']),
           E=0.00025,
           W=0.00025,
           R=0) -> x.F
## Set the prediction distribution equal to the filtering distribution
x.P <- x.F

tstep <- 1
lik <- vector(mode='numeric', length=length(times)-1)
while (tstep < length(times)) {
    ## For each of the Np particles
    for (i in 1:Np) {
        ## generate the food addition events
        events <- subset(eventdat, time >= times[tstep] & time < times[tstep+1])
        events$value <- rnorm(nrow(events),
                              mean=events$value,
                              sd=pars['Ferr'])
        ## obtain a sample of points from the prediction distribution by simulating the model forward
        ode(x.F[i,2:5] %>% unlist,
            times=seq(times[tstep], times[tstep+1]),
            func="derivs",
            parms=pars,
            dllname="debStochEnv",
            initfunc="initmod",
            events=list(data=events)) %>%
            tail(1) -> x.P[i,]
    }
    ## determine the weights by computing the probability of observing the data, given the points in the prediction distribution
    sapply((x.P$W/pars['xi'])^(1/pars['q']),
           function(l)
               dnorm(x=data$length[data$age==times[tstep+1]],
                     mean=l,
                     sd=pars['Lobs'],
                     log=TRUE) %>% sum
           ) +
        sapply(x.P$R,
               function(r)
                   dnorm(x=data$eggs[data$age==times[tstep+1]],
                         mean=r,
                         sd=pars['Robs'],
                         log=TRUE) %>% sum
               ) -> weights

    ## conditional likelihood for this timestep is the mean probability across the points in the prediction distribution
    lik[tstep] <- mean(weights)

    ## use the weights to update the filtering distribution
    sample(order(weights),
           length(weights),
           prob=


    (x[[i]][which(x[[i]][,'T']==t),'W']/pars['xi'])^(1/pars['q'])

sapply(times[2:(tstep+1)], function(t) (x[[i]][t,'W']/pars['xi'])^(1/pars['q']))




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



