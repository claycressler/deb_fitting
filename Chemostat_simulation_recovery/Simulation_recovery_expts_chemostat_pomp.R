## This R-file encodes the DEB models of Cressler_Nelson_style_DEB.R
## into a pomp() object for the purposes of fitting simulated data,
## and also eventually Adriana's data.

require(pomp)

source('DEB_chemostat_constantK_C-code.R')

## Load a simulated dataset
load('~/Dropbox/Cressler/DEB_fitting/LH_chemostat_constantK_fastE_F05.rda')
## Create an "observed" dataset consisting only of length and neonate
## counts at two day intervals with measurement error.
simdata <- out11.fast[seq(11,751,20),c('time','y2','y4')] ## the raw simulated data

## Column names of the dataset should match with what is expected by
## the measurement model for state variable names
colnames(simdata) <- c('Age','Lobs','Robs')
obsdata <- simdata
L_sd <- 0.05 ## length measurement standard deviation
obsdata$Lobs <- rnorm(length(simdata$Lobs),simdata$Lobs,L_sd)
obsdata$Robs <- rpois(length(simdata$Robs),simdata$Robs)

## Build a pomp object
pompBuilder(
            name='DEB_chemostat_constantK',
            data=obsdata,
            times='Age',
            t0=0,
            step.fn=stepfn,
            step.fn.delta.t=0.1,
            dmeasure=dmeas,
            rmeasure=rmeas,
            skeleton.type='vectorfield',
            skeleton=skel,
            parameter.transform=trans,
            parameter.inv.transform=untrans,
            statenames=c('E','L','Re','R'),
            paramnames=c('K','km','eG','eR','v','Rmbar','f','E.0','L.0','Re.0','R.0','L_sd','PA_sd','PC_sd')) -> deb

## Okay, so now I have a compiled C object encoding the DEB model with
## chemostat ingestion and constant kappa. I should be able to specify
## parameters and create simulations using this code directly, rather
## than using the code in
## Simulation_recovery_expt_data_generation.R. I can use those
## parameter values, with one exception: we decided to have a single
## equation for "reproduction energy", rather than two DEs for
## maturity energy and egg production, so the reproduction energy
## equation include eR, which means that Rmbar needs to be modified to
## Rmbar/eR as well.
F <- 0.05
params <- c(K=0.6, km=0.33, eG=0.0017, eR=0.00868, v=18.1, Rmbar=0.0189, f=(0.005187*F/(3.09e-5+F)), E.0=(0.005187/18.1*0.85^3), L.0=0.85, Re.0=0, R.0=0, L_sd=0.05, PA_sd=0, PC_sd=0)
## Add these parameters to the pomp object
deb <- pomp(deb, params=params)
## compute a deterministic trajectory
x <- trajectory(deb, as.data.frame=TRUE, times=seq(0,75,by=0.1))

## Compare the trajectory 'x' to the dataset loaded previously to
## create the first stochastic simulation for trajectory
## matching. They should be, and are, identical.
plot(x$time, x$L, type='l', lwd=2)
points(out11.fast$time, out11.fast$y2, type='l', lty=2, col=gray(1))

## I can now also do some simple trajectory matching.  Let's begin by
## going one parameter at a time, and seeing how well it is able to
## recover single parameters. I will go in order: K, km, eG, eR, v, and Rmbar.
## Modify the initial guess for K away from the true parameter.
## Note that I am fitting data WITH measurement error!
guess <- coef(deb)
guess['K'] <- 0.8 ## truth is 0.6
tm <- traj.match(deb,start=guess,transform=TRUE,est='K',method='subplex',maxit=1000)
coef(tm)['K'] ## 0.598

## Modify km
guess <- coef(deb)
guess['km'] <- 0.6 ## truth is 0.33
tm <- traj.match(deb,start=guess,transform=TRUE,est='km',method='subplex',maxit=1000)
coef(tm)['km'] ## 0.3309

## Modify eG
guess <- coef(deb)
guess['eG'] <- 0.01 ## truth is 0.0017
tm <- traj.match(deb,start=guess,transform=TRUE,est='eG',method='subplex',maxit=1000)
coef(tm)['eG'] ## 0.001704

## Modify eR
guess <- coef(deb)
guess['eR'] <- 0.03 ## truth is 0.00868
tm <- traj.match(deb,start=guess,transform=TRUE,est='eR',method='subplex',maxit=1000)
coef(tm)['eR'] ## 0.008617

## Modify v
guess <- coef(deb)
guess['v'] <- 5 ## truth is 18.1
tm <- traj.match(deb,start=guess,transform=TRUE,est='v',method='subplex',maxit=1000)
coef(tm)['v'] ## 16.23

## Okay, so far so good. What about estimating multiple parameters at
## the same time? K and v, for example?
guess <- coef(deb)
guess[c('K','v')] <- c(0.8, 2) ## truth is 0.6, 18.1
tm <- traj.match(deb,start=guess,transform=TRUE,est=c('K','v'),method='subplex',maxit=1000)
coef(tm)[c('K','v')] ## 0.598, 17.42
## Great - got pretty close, even with the error and the starts far from the correct.

## Try the (nearly) full set of parameters - I will not try to
## estimate the amount of reproduction energy or eggs (set to 0) or
## the length measurement error standard deviation.
guess <- coef(deb)
guess['K'] <- runif(1)
guess['km'] <- runif(1,min=0,max=2)*coef(deb)['km']
guess['eR'] <- runif(1,min=0,max=2)*coef(deb)['eR']
guess['eG'] <- runif(1,min=0,max=2)*coef(deb)['eG']
guess['v'] <- runif(1,min=0,max=2)*coef(deb)['v']
guess['E.0'] <- runif(1,min=0,max=2)*coef(deb)['E.0']
guess['L.0'] <- runif(1,min=0,max=2)*coef(deb)['L.0']
tm <- traj.match(deb,start=guess,transform=TRUE,est=c('K','km','eR','eG','v','E.0','L.0'),method='subplex',maxit=10000)
coef(tm)
## Varying degrees of success here each time the thing is run -
## probably because my initial guess may be in a really uninformative
## part of parameter space. Probably even for trajectory matching it
## would be useful to know where the informative parameter
## combinations are. This means revisiting Aaron and Ed's strategy for
## finding parameter combinations. What if I start relatively nearby?
guess <- coef(deb)
guess[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')] <- 0.9*guess[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')]
tm1 <- traj.match(deb,start=guess,transform=TRUE,est=c('K','km','eR','eG','v','E.0','L.0'),method='subplex',maxit=10000)
round(1-(coef(tm1)[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')]/coef(deb)[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')]),2)

## And here we start a little further below the truth
guess <- coef(deb)
guess[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')] <- 0.8*guess[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')]
tm2 <- traj.match(deb,start=guess,transform=TRUE,est=c('K','km','eR','eG','v','E.0','L.0'),method='subplex',maxit=10000)
round(1-(coef(tm2)[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')]/coef(deb)[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')]),2)

## And here we start a little above the truth
guess <- coef(deb)
guess[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')] <- 1.1*guess[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')]
tm3 <- traj.match(deb,start=guess,transform=TRUE,est=c('K','km','eR','eG','v','E.0','L.0'),method='subplex',maxit=10000)
round(1-(coef(tm3)[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')]/coef(deb)[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')]),2)

## Overall, it does pretty well. The percent that it misses by is
## often less than 10-20%, except for v. Seems likely that v may be
## hard to estimate, given the results of Martin et al. However, I
## started all of these parameter values pretty damn near the
## truth. Likely I will need to use an integrated strategy to get near
## the truth before I begin trajectory matching.

## What if I fit a different dataset? Take the constant kappa with
## chemostat ingestion data for low food and low v value, and see how
## it works there.
load('~/Dropbox/Cressler/DEB_fitting/LH_chemostat_constantK_slowE_F00625.rda')
simdata <- out11.slow[seq(11,751,20),c('time','y2','y4')] ## the raw simulated data
colnames(simdata) <- c('Age','Lobs','Robs')
obsdata <- simdata
L_sd <- 0.05 ## length measurement standard deviation
obsdata$Lobs <- rnorm(length(simdata$Lobs),simdata$Lobs,L_sd)
obsdata$Robs <- rpois(length(simdata$Robs),simdata$Robs)

## Rebuild the pomp object
pompBuilder(
            name='DEB_chemostat_constantK',
            data=obsdata,
            times='Age',
            t0=0,
            step.fn=stepfn,
            step.fn.delta.t=0.1,
            dmeasure=dmeas,
            rmeasure=rmeas,
            skeleton.type='vectorfield',
            skeleton=skel,
            parameter.transform=trans,
            parameter.inv.transform=untrans,
            statenames=c('E','L','Re','R'),
            paramnames=c('K','km','eG','eR','v','Rmbar','f','E.0','L.0','Re.0','R.0','L_sd','PA_sd','PC_sd')) -> deb

## Add the true parameters back in
F <- 0.00625
params <- c(K=0.6, km=0.33, eG=0.0017, eR=0.00868, v=1.81, Rmbar=0.0189, f=(0.005187*F/(3.09e-5+F)), E.0=(0.005187/1.81*0.85^3), L.0=0.85, Re.0=0, R.0=0, L_sd=0.05, PA_sd=0, PC_sd=0)
deb <- pomp(deb, params=params)

## Compute the deterministic trajectory and compare - again, they are
## identical
x <- trajectory(deb, as.data.frame=T, times=seq(0,75,0.1))
plot(x$time, x$L, type='l', lwd=2)
points(out11.slow$time, out11.slow$y2, type='l', lty=2, col=gray(1))

## Use trajectory matching as before, parameter guesses just below the truth
guess <- coef(deb)
guess[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')] <- 0.8*guess[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')]
tm4 <- traj.match(deb,start=guess,transform=TRUE,est=c('K','km','eR','eG','v','E.0','L.0'),method='subplex',maxit=10000)
round(1-(coef(tm4)[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')]/coef(deb)[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')]),3)

## Use trajectory matching as before, parameter guesses just above the truth
guess <- coef(deb)
guess[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')] <- 1.2*guess[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')]
tm5 <- traj.match(deb,start=guess,transform=TRUE,est=c('K','km','eR','eG','v','E.0','L.0'),method='subplex',maxit=10000)
round(1-(coef(tm5)[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')]/coef(deb)[c('K','km','eR','eG','v','Rmbar','f','E.0','L.0')]),3)

## Again, not too shabby! I can get fairly close, provided that I
## start in a reasonable neighborhood of the truth. The parameters
## that are missing by wider margins are those that are teeny-tiny in
## absolute value (eR, eG, f, and E.0, in particular), so these large
## proportional differences probably don't matter much for the actual
## dynamics.



## Create my own - in this case, add food every two days
derivfun <- function (t,y,parms)
    list (-0.05 * y)

rootfun <- function (t,y,parms)
    return(sin(pi/2*t))

eventfun <- function(t,y,parms)
    return(y + runif(1))

yini <- 0.8
times <- seq(0,20,0.1)

out <- lsodar(func=derivfun, y = yini, times=times,
              rootfunc = rootfun, events = list(func=eventfun, root = TRUE))

plot(out, type = "l", lwd = 2, main = "lsodar with event")


## I can now also simulate stochastic trajectories of this model.
## compute a stochastic realization - since PA_sd=PC_sd=0, this should
## be identical to the deterministic trajectory.
#y <- simulate(deb, nsim=1)

