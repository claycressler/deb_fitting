source("Growth_reproduction_trajectory_fitting_functions_2.R")

pars <- c(Imax=22500, fh=10000, g=1.45, rho=0.1, eps=44.5e-9, V=30, F0=1000000/30, xi=2.62e-3, q=2.4, K=0.3, km=0.15, ER=1.51e-3, v=10, Lobs=0.1)
y0 <- c(F=1000000/30, E=0.00025, W=0.00025, R=0)
times <- seq(0, 35, 0.001)
## feeding schedule
eventdat <- data.frame(var="F",
                       time=1:35,
                       value=unname(pars["F0"]),
                       method=rep(c(rep("add",4),"rep"),max(times)/5))
ode(y0, times=0:35, func="derivs", parms=pars, dllname="deb2", initfunc="initmod", events=list(data=eventdat)) %>% as.data.frame -> out2
days <- c(5,10,12,15,18,25,30,35)
set.seed(1239478)
data.frame(times=rep(days, each=12),
           length=sapply(with(out2, ((E[days+1]+W[days+1])/pars["xi"])^(1/pars["q"])),
               function(x)
                   rnorm(12, mean=x, sd=pars["Lobs"])
                         ) %>% as.numeric,
           eggs=sapply(with(out2, R-R[which((E+W) < 0.005) %>% max])[days+1],
               function(x)
                   rpois(12, lambda=x)
                       ) %>% as.numeric
           ) -> data
data$eggs[is.na(data$eggs)] <- 0
saveRDS(data, file="simulated_dataset.RDS")

parorder <- c("Imax","fh","g","rho","eps","V","F0","xi","q","K","km","ER","v","Lobs")
fixpars <- pars[c("Imax","fh","g","rho","eps","V","F0","xi","q")]
estpars <- pars[c("K","km","ER","v","Lobs")]
transform <- c("logit", rep("log",4))

## testing
obj(estpars=par_transform(estpars, transform),
    data=data,
    fixpars=fixpars,
    parorder=parorder,
    transform=transform)
traj_match(estpars, fixpars, parorder, transform, data, eval.only=TRUE)

## can we recover all of the parameters if we treat fh as known?
box <- cbind(lower=c(K=0, km=0.001, ER=0.00001, v=0.1, Lobs=0.0001),
             upper=c(K=1, km=10, ER=1, v=100, Lobs=2))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=250000) %>%
    apply(., 1, as.list) %>%
        lapply(., unlist) -> guesses

t1 <- Sys.time()
mclapply(guesses,
         traj_match,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         eval.only=TRUE,
         mc.cores=10) %>%
    lapply(., function(x) x$lik) %>%
        unlist -> guess_lik
t2 <- Sys.time()
print(t2-t1)
guesses[order(guess_lik)[1:1000]] -> refine
mclapply(refine,
         traj_match,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         eval.only=FALSE,
         method="Nelder-Mead",
         mc.cores=10) -> refine_lik
refine_lik %>%
    lapply(., unlist) %>%
        unlist %>%
            matrix(., ncol=(nrow(box)+2), byrow=TRUE, dimnames=list(1:length(refine_lik), c(rownames(box),"lik","conv"))) %>%
                as.data.frame -> refine_pars
saveRDS(refine_pars, file="Growth_reproduction_trajectory_fitting_dyn_food.RDS")

####################################################
## Do the estimates tighten up if you fix v and fit the other
## parameters, as we saw before?
source("Growth_reproduction_trajectory_fitting_functions_2.R")
data <- readRDS("simulated_dataset.RDS")
transform <- c("logit", rep("log",3))
parorder <- c("Imax","fh","g","rho","eps","V","F0","xi","q","K","km","ER","v","Lobs")

## initial guesses for the estimated parameters
box <- cbind(lower=c(K=0, km=0.001, ER=0.00001, Lobs=0.0001),
             upper=c(K=1, km=10, ER=1, Lobs=2))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=250000) %>%
    apply(., 1, as.list) %>%
        lapply(., unlist) -> guesses

v_vals <- c(seq(10, 100, 10), seq(200, 1000, 100))
ests <- vector(mode='list', length=length(v_vals))
for (i in 1:length(ests)) {
    print(i)
    fixpars <- c(Imax=22500, fh=10000, g=1.45, rho=0.1, eps=44.5e-9, V=30, F0=1000000/30, xi=2.62e-3, q=2.4, v=v_vals[i])
    eventdat <- data.frame(var="F",
                           time=1:35,
                           value=unname(fixpars["F0"]),
                           method=rep(c(rep("add",4),"rep"),35/5))
    t1 <- Sys.time()
    mclapply(guesses,
             traj_match,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=data,
             eval.only=TRUE,
             mc.cores=10) %>%
        lapply(., function(x) x$lik) %>%
            unlist -> guess_lik
    t2 <- Sys.time()
    print(t2-t1)

    t1 <- Sys.time()
    guesses[order(guess_lik)[1:1000]] -> refine
    mclapply(refine,
             traj_match,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=data,
             eval.only=FALSE,
             method="Nelder-Mead",
             mc.cores=10) -> refine_lik
    t2 <- Sys.time()
    print(t2-t1)

    refine_lik %>%
        lapply(., unlist) %>%
            unlist %>%
                matrix(., ncol=(nrow(box)+2), byrow=TRUE, dimnames=list(1:length(refine_lik), c(rownames(box),"lik","conv"))) %>%
                    as.data.frame -> refine_pars
    ests[[i]] <- refine_pars

    saveRDS(ests, file="Growth_reproduction_trajectory_fitting_dyn_food_profile_v.RDS")
}

###################################################
### Clearly, v is still not identifiable, so simply pick a value for it and use that. However, what if you make the inference problem harder by trying to estimate a few more parameters? In particular, what if you need to estimate fh?

x <- readRDS("Feeding_model_fitting_uninfecteds.RDS")
## extract best fitting parameter estimate for every fh value
lapply(x, function(i) subset(i, lik==min(lik))) %>%
    unlist %>%
        matrix(., ncol=6, byrow=TRUE, dimnames=list(1:length(x), colnames(x[[1]]))) %>% as.data.frame -> z
z$fh <- seq(100, 50000, 100)
## linear relationship between fh and Imax that makes it easy to set Imax for
## the challenge is that Imax and g vary with the value of fh.
## Imax has a very simple relationship with fh - a perfectly linear function
## lm(Imax~fh) has a R^2 of 1, with coefs
with(z, summary(lm(Imax~fh)))
calc_Imax <- function(fh) 10616.500976 + 1.129421*fh

## fh and g have a nonlinear relationship, but a cubic function predicts it with R^2 of 0.99
with(z, lm(g~fh+I(fh^2)+I(fh^3)) %>% summary)
calc_g <- function(fh) 1.38413 + 9.839020e-6*fh - 2.738144e-10*fh^2 + 2.648817e-15*fh^3


###### Simulation-recovery exercise, again
source("Growth_reproduction_trajectory_fitting_functions_3.R")
data <- readRDS("simulated_dataset.RDS")

## initial guesses for the estimated parameters
box <- cbind(lower=c(fh=100, K=0, km=0.001, ER=0.00001, Lobs=0.0001),
             upper=c(fh=50000, K=1, km=10, ER=1, Lobs=2))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=250000) %>%
    apply(., 1, as.list) %>%
        lapply(., unlist) -> guesses
transform <- c("log", "logit", rep("log",3))
parorder <- c("Imax","fh","g","rho","eps","V","F0","xi","q","K","km","ER","v","Lobs")

v_vals <- c(seq(10, 100, 10), seq(200, 1000, 100))
ests <- vector(mode='list', length=length(v_vals))
for (i in 1:length(ests)) {
    print(i)
    ## Imax and g are included in fixpars even though they will actually be
    ## set by the estimate of
    fixpars <- c(Imax=22500, g=1.45, rho=0.1, eps=44.5e-9, V=30, F0=1000000/30, xi=2.62e-3, q=2.4, v=v_vals[i])
    t1 <- Sys.time()
    eventdat <- data.frame(var="F",
                           time=1:35,
                           value=unname(fixpars["F0"]),
                           method=rep(c(rep("add",4),"rep"),35/5))
    mclapply(guesses,
             traj_match,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=data,
             events=eventdat,
             eval.only=TRUE,
             mc.cores=10) %>%
        lapply(., function(x) x$lik) %>%
            unlist -> guess_lik
    t2 <- Sys.time()
    print(t2-t1)

    t1 <- Sys.time()
    guesses[order(guess_lik)[1:1000]] -> refine
    mclapply(refine,
             traj_match,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=data,
             events=eventdat,
             eval.only=FALSE,
             method="Nelder-Mead",
             mc.cores=10) -> refine_lik
    t2 <- Sys.time()
    print(t2-t1)

    refine_lik %>%
        lapply(., unlist) %>%
            unlist %>%
                matrix(., ncol=(nrow(box)+2), byrow=TRUE, dimnames=list(1:length(refine_lik), c(rownames(box),"lik","conv"))) %>%
                    as.data.frame -> refine_pars
    ests[[i]] <- refine_pars

    saveRDS(ests, file="~/Dropbox/Growth_reproduction_trajectory_fitting_dyn_food_profile_v_2.RDS")
}

########################################
## What about estimating the assimilation efficiency? Can you do that while simultaneously estimating fh?
calc_Imax <- function(fh) 10616.500976 + 1.129421*fh
calc_g <- function(fh) 1.38413 + 9.839020e-6*fh - 2.738144e-10*fh^2 + 2.648817e-15*fh^3

source("Growth_reproduction_trajectory_fitting_functions_3.R")
data <- readRDS("simulated_dataset.RDS")

## initial guesses for the estimated parameters
box <- cbind(lower=c(fh=100, rho=0, K=0, km=0.001, ER=0.00001, Lobs=0.0001),
             upper=c(fh=50000, rho=1, K=1, km=10, ER=1, Lobs=2))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=250000) %>%
    apply(., 1, as.list) %>%
        lapply(., unlist) -> guesses
transform <- c("log", rep("logit",2), rep("log",3))
parorder <- c("Imax","fh","g","rho","eps","V","F0","xi","q","K","km","ER","v","Lobs")

fixpars <- c(Imax=22500, g=1.45, eps=44.5e-9, V=30, F0=1000000/30, xi=2.62e-3, q=2.4, v=100)
eventdat <- data.frame(var="F",
                       time=1:35,
                       value=unname(fixpars["F0"]),
                       method=rep(c(rep("add",4),"rep"),35/5))
mclapply(guesses,
         traj_match,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         events=eventdat,
         eval.only=TRUE,
         mc.cores=10) %>%
    lapply(., function(x) x$lik) %>%
        unlist -> guess_lik
guesses[order(guess_lik)[1:1000]] -> refine
mclapply(refine,
         traj_match,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         events=eventdat,
         eval.only=FALSE,
         method="Nelder-Mead",
         mc.cores=10) -> refine_lik
refine_lik %>%
    lapply(., unlist) %>%
        unlist %>%
            matrix(., ncol=(nrow(box)+2), byrow=TRUE, dimnames=list(1:length(refine_lik), c(rownames(box),"lik","conv"))) %>%
                as.data.frame -> refine_pars
saveRDS(refine_pars, file="Estimating-rho.RDS")

############################################
##### Profile likelihood for rho
source("Growth_reproduction_trajectory_fitting_functions_3.R")
data <- readRDS("simulated_dataset.RDS")

## initial guesses for the estimated parameters
box <- cbind(lower=c(fh=100, K=0, km=0.001, ER=0.00001, Lobs=0.001),
             upper=c(fh=50000, K=1, km=1, ER=0.1, Lobs=2))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=250000) %>%
    apply(., 1, as.list) %>%
        lapply(., unlist) -> guesses
transform <- c("log", "logit", rep("log",3))
parorder <- c("Imax","fh","g","rho","eps","V","F0","xi","q","K","km","ER","v","Lobs")



rho_vals <- seq(0.1, 0.9, 0.05)
ests <- vector(mode='list', length=length(rho_vals))
for (i in 10:length(ests)) {
    print(i)
    ## Imax and g are included in fixpars even though they will actually be
    ## set by the estimate of
    fixpars <- c(Imax=22500, g=1.45, rho=rho_vals[i], eps=44.5e-9, V=30, F0=1000000/30, xi=2.62e-3, q=2.4, v=100)
    eventdat <- data.frame(var="F",
                           time=1:35,
                           value=unname(fixpars["F0"]),
                           method=rep(c(rep("add",4),"rep"),35/5))
    t1 <- Sys.time()
    mclapply(guesses,
             traj_match,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=data,
             events=eventdat,
             eval.only=TRUE,
             mc.cores=15) %>%
        lapply(., function(x) x$lik) %>%
            unlist -> guess_lik
    t2 <- Sys.time()
    print(t2-t1)

    t1 <- Sys.time()
    guesses[order(guess_lik)[1:500]] -> refine
    mclapply(refine,
             traj_match,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=data,
             events=eventdat,
             eval.only=FALSE,
             method="Nelder-Mead",
             mc.cores=15) -> refine_lik
    t2 <- Sys.time()
    print(t2-t1)

    refine_lik %>%
        lapply(., unlist) %>%
            unlist %>%
                matrix(., ncol=(nrow(box)+2), byrow=TRUE, dimnames=list(1:length(refine_lik), c(rownames(box),"lik","conv"))) %>%
                    as.data.frame -> refine_pars
    ests[[i]] <- refine_pars

    saveRDS(ests, file="~/Dropbox/Growth_reproduction_trajectory_fitting_dyn_food_profile_rho.RDS")
}

############################################################

source("Growth_reproduction_trajectory_fitting_functions_3.R")
data <- readRDS("simulated_dataset.RDS")

## initial guesses for the estimated parameters
box <- cbind(lower=c(fh=100, rho=0, K=0, km=0.001, Lobs=0.0001),
             upper=c(fh=50000, rho=1, K=1, km=10, Lobs=2))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=250000) %>%
    apply(., 1, as.list) %>%
        lapply(., unlist) -> guesses
transform <- c("log", rep("logit",2), rep("log",2))
parorder <- c("Imax","fh","g","rho","eps","V","F0","xi","q","K","km","ER","v","Lobs")

v_vals <- c(seq(10, 100, 10), seq(200, 1000, 100))
ests <- vector(mode='list', length=length(v_vals))
for (i in 6:(length(ests)-1)) {
    print(i)
    fixpars <- c(Imax=22500, g=1.45, eps=44.5e-9, V=30, F0=1000000/30, xi=2.62e-3, q=2.4, ER=1.51e-3, v=v_vals[i])
    eventdat <- data.frame(var="F",
                           time=1:35,
                           value=unname(fixpars["F0"]),
                           method=rep(c(rep("add",4),"rep"),35/5))
    mclapply(guesses,
             traj_match,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=data,
             events=eventdat,
             eval.only=TRUE,
             mc.cores=12) %>%
        lapply(., function(x) x$lik) %>%
            unlist -> guess_lik
    guesses[order(guess_lik)[1:500]] -> refine
    mclapply(refine,
             traj_match,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=data,
             events=eventdat,
             eval.only=FALSE,
             method="Nelder-Mead",
             mc.cores=12) -> refine_lik
    refine_lik %>%
        lapply(., unlist) %>%
            unlist %>%
                matrix(., ncol=(nrow(box)+2), byrow=TRUE, dimnames=list(1:length(refine_lik), c(rownames(box),"lik","conv"))) %>%
                    as.data.frame -> refine_pars
    ests[[i]] <- refine_pars
    saveRDS(ests, file="Growth_reproduction_trajectory_fitting_dyn_food_profile_v_est_rho.RDS")
}


############################################################
## Generate some random parameter sets and observed datasets to make sure the conclusions from the above hold generally. Specifically, it appears that it is very possible to estimate the values of fh, rho, K, km, and Lobs, given a fixed value of ER and v.

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
box <- cbind(lower=c(fh=100, rho=0, K=0, km=0.001, Lobs=0.0001),
             upper=c(fh=50000, rho=1, K=1, km=10, Lobs=2))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=250000) %>%
    apply(., 1, as.list) %>%
        lapply(., unlist) -> guesses
transform <- c("log", rep("logit",2), rep("log",2))
parorder <- c("Imax","fh","g","rho","eps","V","F0","xi","q","K","km","ER","v","Lobs")

ests <- vector(mode='list', length=length(datasets))
for (i in 1:length(datasets)) {
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
    saveRDS(ests, file="~/Dropbox/Growth_reproduction_trajectory_fitting_dyn_food_multiple_datasets.RDS")
}

ests <- vector(mode='list', length=length(datasets))
for (i in 1:length(datasets)) {
    print(i)
    ## although ER varied in the simulations, we need to fix its value
    ## for the recovery to have any hope of estimating rho (and
    ## kappa). Similarly, we are holding v fixed at 100, since it
    ## doesn't seem to possible to recover this parameter and it
    ## doesn't seem to really affect any of the other parameter
    ## estimatesa
    fixpars <- c(Imax=22500, g=1.45, eps=44.5e-9, V=30, F0=1000000/30, xi=2.62e-3, q=2.4, datasets[[i]]$params["ER"], v=10)
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
    saveRDS(ests, file="~/Dropbox/Growth_reproduction_trajectory_fitting_dyn_food_multiple_datasets_v=10.RDS")
}
