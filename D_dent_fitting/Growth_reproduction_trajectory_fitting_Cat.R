source("Growth_reproduction_trajectory_fitting_functions_Cat.R")

## Load Cat's data
x <- read.csv("Cat_data/uninfected_growth_reproduction.csv")
data <- x[1:103,1:3]

pars <- c(Imax=calc_Imax(10000), # fixed based on feeding data fitting and fh value
          fh=10000, # estimated from fitting growth/reproduction data
          g=calc_g(10000), # fixed based on feeding data fitting and fh value
          rho=0.2, # estimated from fitting growth/reproduction data
          eps=44.5e-9, # fixed based on measured carbon content of algae
          V=30, # fixed based on experimental conditions
          F0=1000000/30, # fixed based on experimental conditions
          xi=1.8e-3, # fixed based on Hall et al. 2009 length-weight regression
          q=3, # fixed based on Hall et al. 2009 length-weight regression
          K=0.5, # estimated from growth/reproduction data
          km=0.1, # estimated from growth/reproduction data
          ER=1.51e-3, # fixed based on observed mass of a neonate
          v=10, # fixed based on the fact that it doesn't affect the fitting
          Lobs=0.1,
          Winit=0.00225,
          Wmat=0.005) # estimated from growth/reproduction data

## Feeding schedule
eventdat <- data.frame(var="F",
                       time=1:35,
                       value=unname(pars["F0"]),
                       method=rep(c(rep("add",4),"rep"),35/5))
## Initial conditions
y0 <- c(F=unname(pars["F0"]), E=unname(pars["Winit"])/2, W=unname(pars["Winit"])/2, R=0)
## Simulate at these parameters
ode(y0, times=0:35, func="derivs", parms=pars, dllname="deb_Cat", initfunc="initmod", events=list(data=eventdat)) %>% as.data.frame -> out
## Subtract reproduction prior to sexual maturity and calculate observed length
mutate(out,
       L=((E+W)/pars["xi"])^(1/pars["q"]),
       R=R-R[which((E+W) < pars["Wmat"]) %>% max]
       ) -> out
out$R[out$R < 0] <- 0

parorder <- c("Imax","fh","g","rho","eps","V","F0","xi","q","K","km","ER","v","Lobs","Winit","Wmat")
fixpars <- pars[c("Imax","g","eps","V","F0","xi","q","ER","v","Winit","Wmat")]
estpars <- pars[c("fh","rho","K","km","Lobs")]
transform <- c("log",rep("logit",2), rep("log",2))

## testing
obj(estpars=par_transform(estpars, transform),
    data=data,
    fixpars=fixpars,
    parorder=parorder,
    transform=transform,
    events=eventdat)
traj_match(estpars, fixpars, parorder, transform, data, events=eventdat, eval.only=TRUE)
## Calculate the probability of seeing the real data, if these were
## the true parameter values
(dnorm(data$length,
       sapply(data$times,
              function(t) out$L[out$time==t]),
       unname(pars["Lobs"]),
       log=TRUE)  +
     dpois(data$eggs,
           sapply(data$times,
                  function(t) out$R[out$time==t]),
           log=TRUE)
 ) %>% sum(., na.rm=TRUE) -> lik

#########################################################################
#########################################################################
###### FITTING GROWTH AND REPRODUCTION DATA FOR UNINFECTED ANIMALS ######
#########################################################################
#########################################################################

parorder <- c("Imax","fh","g","rho","eps","V","F0","xi","q","K","km","ER","v","Lobs","Winit","Wmat")
transform <- c("log",rep("logit",2), rep("log",2))
pars <- c(Imax=calc_Imax(10000), # fixed based on feeding data fitting and fh value
          fh=10000, # estimated from fitting growth/reproduction data
          g=calc_g(10000), # fixed based on feeding data fitting and fh value
          rho=0.2, # estimated from fitting growth/reproduction data
          eps=44.5e-9, # fixed based on measured carbon content of algae
          V=30, # fixed based on experimental conditions
          F0=1000000/30, # fixed based on experimental conditions
          xi=1.8e-3, # fixed based on Hall et al. 2009 length-weight regression
          q=3, # fixed based on Hall et al. 2009 length-weight regression
          K=0.5, # estimated from growth/reproduction data
          km=0.1, # estimated from growth/reproduction data
          ER=1.51e-3, # fixed based on observed mass of a neonate
          v=10, # fixed based on the fact that it doesn't affect the fitting
          Lobs=0.1,
          Winit=0.00225,
          Wmat=0.005) # estimated from growth/reproduction data

## initial guesses for the estimated parameters
box <- cbind(lower=c(fh=100, rho=0, K=0, km=0.001, Lobs=0.0001),
             upper=c(fh=50000, rho=1, K=1, km=10, Lobs=2))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=250000) %>%
    apply(., 1, as.list) %>%
        lapply(., unlist) -> guesses

ests <- vector(mode='list', length=length(datasets))
for (i in 1:length(datasets)) {
    print(i)
    ## although ER varied in the simulations, we need to fix its value
    ## for the recovery to have any hope of estimating rho (and
    ## kappa). Similarly, we are holding v fixed at 100, since it
    ## doesn't seem to possible to recover this parameter and it
    ## doesn't seem to really affect any of the other parameter
    ## estimatesa
    fixpars <- pars[c("Imax","g","eps","V","F0","xi","q","ER","v","Winit","Wmat")]
    eventdat <- data.frame(var="F",
                           time=1:35,
                           value=unname(fixpars["F0"]),
                           method=rep(c(rep("add",4),"rep"),35/5))
    mclapply(guesses[1:20],
             traj_match,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=data,
             events=eventdat,
             eval.only=TRUE,
             mc.cores=2) %>%
        lapply(., function(x) x$lik) %>%
            unlist -> guess_lik
    guesses[order(guess_lik)[1:1000]] -> refine
    mclapply(refine,
             traj_match,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=datasets[[i]]$data,
             events=eventdat,
             eval.only=FALSE,
             method="Nelder-Mead",
             mc.cores=12) -> refine_lik
    refine_lik %>%
        lapply(., unlist) %>%
            unlist %>%
                matrix(., ncol=(nrow(box)+2), byrow=TRUE, dimnames=list(1:length(refine_lik), c(rownames(box),"lik","conv"))) %>%
                    as.data.frame -> refine_pars
    refine_pars <- arrange(refine_pars, lik)
    print(refine_pars[1,1:5]-datasets[[i]]$params[c("fh","rho","K","km","Lobs")])

    ests[[i]] <- refine_pars
    saveRDS(ests, file="~/Dropbox/Growth_reproduction_trajectory_fitting_dyn_food_multiple_datasets_take_2.RDS")
}

##############################################################
###### Refitting some of the datasets with more constrained
###### initial guess of the parameters
##############################################################

source("Growth_reproduction_trajectory_fitting_functions_3.R")
pars <- c(Imax=22500, fh=10000, g=1.45, rho=0.5, eps=44.5e-9, V=30, F0=1000000/30, xi=2.62e-3, q=2.4, K=0.5, km=0.3, ER=1.51e-3, v=10, Lobs=0.1)
parorder <- c("Imax","fh","g","rho","eps","V","F0","xi","q","K","km","ER","v","Lobs")
fixpars <- pars[c("Imax","g","eps","V","F0","xi","q")]
varpars <- pars[c("fh","rho","K","km","ER","v","Lobs")]
y0 <- c(F=1000000/30, E=0.00025, W=0.00025, R=0)
times <- seq(0, 35, 0.001)
## feeding schedule
eventdat <- data.frame(var="F",
                       time=1:35,
                       value=unname(pars["F0"]),
                       method=rep(c(rep("add",4),"rep"),max(times)/5))
days <- c(5,10,12,15,18,25,30,35)
datasets <- vector(mode='list', length=20)
set.seed(1231239478)
for (i in 1:20) {
    data.frame(times=rep(days, each=12),
               length=rep(0,96),
               eggs=rep(0,96)) -> data
    ## ensure that the data give you something that looks reasonably like a Daphnia
    while ( (min(data$length) < 0.5) | (max(data$length) > 5) |
               any(subset(data, times==5)$length > 2) |
                   coef(lm(length~times,data))[2] < 0.01 |
                       summary(lm(length~times, data))$coefficients[2,4] > 0.05  |
                           (max(data[,"eggs"]) < 5) | (max(data[,"eggs"]) > 200) |
                               any(subset(data, times==5)$eggs > 0) |
                                   any(is.na(data))
           ) {
        ## GENERATE NOVEL PARAMETERS
        rnorm(length(varpars), mean=varpars, sd=varpars/2) -> p
        names(p) <- names(varpars)
        while (p["K"] > 0.9 | p["K"] < 0.1 | any(p < 0) | p["rho"] > 0.9 | p["rho"] < 0.1) {
            rnorm(length(varpars), mean=varpars, sd=varpars/2) -> p
            names(p) <- names(varpars)
        }
        ## combine fixpars with p
        pars <- c(p, fixpars)
        pars <- pars[match(parorder, names(pars))]
        ## calculate the value of Imax and g, given the value of fh
        pars["Imax"] <- calc_Imax(unname(pars["fh"]))
        pars["g"] <- calc_g(unname(pars["fh"]))
        ## simulate
        ode(y0, times=0:35, func="derivs", parms=pars, dllname="deb2", initfunc="initmod", events=list(data=eventdat)) %>% as.data.frame -> out
        mutate(out, R=R-R[which((E+W) < 5e-3) %>% max]) -> out
        out$R[out$R < 0] <- 0
        data.frame(times=rep(days, each=12),
                   length=sapply(with(out, ((E[days+1]+W[days+1])/fixpars["xi"])^(1/fixpars["q"])),
                       function(x)
                           rnorm(12, mean=x, sd=varpars["Lobs"])
                                 ) %>% as.numeric,
                   eggs=sapply(out$R[days+1],
                       function(x)
                           rpois(12, lambda=x)
                               ) %>% as.numeric
                   ) -> data
    }
    datasets[[i]] <- list(params=pars, data=data)
}

## initial guesses for the estimated parameters
box <- cbind(lower=c(fh=2000, rho=0, K=0, km=0.001, Lobs=0.001),
             upper=c(fh=20000, rho=1, K=1, km=1, Lobs=2))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=250000) %>%
    apply(., 1, as.list) %>%
        lapply(., unlist) -> guesses
transform <- c("log", rep("logit",2), rep("log",2))
parorder <- c("Imax","fh","g","rho","eps","V","F0","xi","q","K","km","ER","v","Lobs")

ests <- vector(mode='list', length=length(datasets))
for (i in c(3,8,9,10,12,17,18,20)) {
    print(i)
    ## although ER varied in the simulations, we need to fix its value
    ## for the recovery to have any hope of estimating rho (and
    ## kappa). Similarly, we are holding v fixed at 100, since it
    ## doesn't seem to possible to recover this parameter and it
    ## doesn't seem to really affect any of the other parameter
    ## estimatesa
    fixpars <- c(Imax=22500, g=1.45, eps=44.5e-9, V=30, F0=1000000/30, xi=2.62e-3, q=2.4, datasets[[i]]$params["ER"], v=100)
    eventdat <- data.frame(var="F",
                           time=1:35,
                           value=unname(fixpars["F0"]),
                           method=rep(c(rep("add",4),"rep"),35/5))
    mclapply(guesses,
             traj_match,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=datasets[[i]]$data,
             events=eventdat,
             eval.only=TRUE,
             mc.cores=12) %>%
        lapply(., function(x) x$lik) %>%
            unlist -> guess_lik
    guesses[order(guess_lik)[1:1000]] -> refine
    mclapply(refine,
             traj_match,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=datasets[[i]]$data,
             events=eventdat,
             eval.only=FALSE,
             method="Nelder-Mead",
             mc.cores=12) -> refine_lik
    refine_lik %>%
        lapply(., unlist) %>%
            unlist %>%
                matrix(., ncol=(nrow(box)+2), byrow=TRUE, dimnames=list(1:length(refine_lik), c(rownames(box),"lik","conv"))) %>%
                    as.data.frame -> refine_pars
    refine_pars <- arrange(refine_pars, lik)
    print(refine_pars[1,1:5]-datasets[[i]]$params[c("fh","rho","K","km","Lobs")])

    ests[[i]] <- refine_pars
    saveRDS(ests, file="~/Dropbox/Growth_reproduction_trajectory_fitting_dyn_food_multiple_datasets_2.RDS")
}
