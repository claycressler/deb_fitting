## Load the code for performing trajectory matching
source("Growth_reproduction_trajectory_matching_real_data_2.R")

## Load Cat's data
x <- read.csv("Cat_data/uninfected_growth_reproduction.csv")
data <- x[1:103,1:3]
data$eggs[is.na(data$eggs)] <- 0 ## set reproduction = 0 for all individuals that have not yet matured
## change the name of 'times' to 'age'
colnames(data)[1] <- 'age'

## I need to think carefully about how to deal with maturation. The DEB model will allow individuals to reproduce immediately. I could treat size at maturity as a parameter to be estimated. Or I could simply say, for each individual, that size at maturity is whatever size reproduction was first observed. I will have to play around with both of these models. I will start with a model where age at maturity has to be estimated along with the other parameters.

## Initial guesses at parameter values
fixpars <- c(Imax=22500, g=1.45, v=10, F0=1e6/30)
estpars <- c(Fh=10000, # estimated from fitting growth/reproduction data
             rho=0.2, # estimated from fitting growth/reproduction data
             K=0.5, # estimated from growth/reproduction data
             km=0.1, # estimated from growth/reproduction data
             ER=1.51e-3, # fixed based on observed mass of a neonate
             v=10, # fixed based on the fact that it doesn't affect the fitting
             Lobs=0.1,
             Robs=2,
             Wmat=0.005) # estimated from growth/reproduction data
transform <- c("log", rep("logit",2), rep("log",6))
parorder <- c("Imax","Fh","g","rho","K","km","v","ER","F0","Lobs","Robs","Wmat")

fixpars["Imax"] <- calc_Imax(unname(estpars["Fh"]))
fixpars["g"] <- calc_g(unname(estpars["Fh"]))
## combine the parameters to be estimated and the fixed parameters
## into a single vector, with an order specified by parorder
pars <- c(estpars, fixpars)
pars[match(parorder, names(pars))] -> pars

## Simulate the model to make sure everything is working
## Set the initial conditions
y0 <- c(F=unname(pars["F0"]),
        E=0,
        W=unname(pars["ER"]/(1+pars["rho"]/pars["v"])),
        M=0,
        R=0)
y0["E"] <- unname(y0["W"]*pars["rho"]/pars["v"])

## feeding schedule - amount of food added each time is stochastic
set.seed(101)
eventdat <- data.frame(var="F",
                       time=1:35,
                       value=unname(pars["F0"]),
                       method=rep(c(rep("add",4),"rep"),7))
## simulate the system
ode(y0,
    times=seq(0,35,0.1),
    func="derivs",
    parms=pars,
    dllname="tm_deb_2",
    initfunc="initmod",
    events=list(data=eventdat)) -> out

## compute the likelihood
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
             dnorm(x=data$eggs[data$age==d],
                   mean=pred$R[pred$time==d],
                   sd=pars["Robs"],
                   log=TRUE) %>% sum
             ) %>% sum
       ) %>% sum -> lik

## compute the likelihood using the 'tm_obj' function
tm_obj(estpars=par_transform(estpars, transform), data=data, fixpars=fixpars, parorder=parorder, transform=transform)

## compute the likelihood using the 'optimizer' function
optimizer(estpars=estpars, fixpars=fixpars, parorder=parorder, transform=transform, obsdata=data, eval.only=TRUE)


####################################################################
####################################################################
########    TRAJECTORY MATCHING OF THE REAL DATA!!!!    ############
####################################################################
####################################################################

## Load the code for performing trajectory matching
source("Growth_reproduction_trajectory_matching_real_data_2.R")

## Load Cat's data
x <- read.csv("Cat_data/uninfected_growth_reproduction.csv")
data <- x[1:103,1:3]
data$eggs[is.na(data$eggs)] <- 0 ## set reproduction = 0 for all individuals that have not yet matured
## change the name of 'times' to 'age'
colnames(data)[1] <- 'age'

## Begin the funnel of opimization by estimating the likelihood of a huge number of parameter guesses
box <- cbind(lower=c(Fh=2000, rho=0, K=0, km=0.001, ER=0.00001, Lobs=0.001, Robs=0.01, Wmat=0.0001),
             upper=c(Fh=20000, rho=1, K=1, km=1, ER=0.01, Lobs=2, Robs=20, Wmat=0.01))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=300000) %>%
                apply(., 1, as.list) %>%
                    lapply(., unlist) -> guesses

## Fixed parameter values (actually, Imax and g are set by the value of Fh)
fixpars <- c(Imax=22500, g=1.45, v=10, F0=1e6/30)
## Parameter transformation to the unconstrained scale
transform <- c("log", rep("logit",2), rep("log",5))
## Order that tm_deb.c expects the parameters to be in
parorder <- c("Imax","Fh","g","rho","K","km","v","ER","F0","Lobs","Robs","Wmat")

## Estimate likelihood
mclapply(guesses,
         optimizer,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         eval.only=TRUE,
         type="trajectory_matching",
         mc.cores=15) %>%
             lapply(., function(x) x$lik) %>%
                 unlist -> guess_lik
guesses[order(guess_lik)[1:3000]] -> refine
mclapply(refine,
         optimizer,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         eval.only=FALSE,
         type='trajectory_matching',
         method='Nelder-Mead',
         mc.cores=15) -> refine_lik
refine_lik %>%
    lapply(., unlist) %>%
        unlist %>%
            matrix(., ncol=(nrow(box)+2), byrow=TRUE, dimnames=list(1:length(refine_lik), c(rownames(box),"lik","conv"))) %>%
                as.data.frame -> refine_pars
refine_pars[order(refine_pars$lik),] -> refine_pars
saveRDS(refine_pars, file="Trajectory_matching_estimates_uninfected_animals_11-29.RDS")
