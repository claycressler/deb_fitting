source("Infected_feeding_model_fitting_functions.R")

data <- read.csv("Cat_data/feeding_infection.csv")
sdata <- subset(data, infected==1 & !is.na(Z))

## extract some parameters from previous fits to remove those from the estimation procedure
readRDS("Feeding_model_fitting_uninfecteds.RDS") %>%
    lapply(., function(l) subset(l, lik==min(lik))[c("F0","Fobs")]) %>%
        unlist %>% matrix(., ncol=2, byrow=T) %>% apply(., 2, mean) -> p
fixpars <- c(fh=1000, F0=p[1], Fobs=p[2])
estpars <- c(a=1, h=2, Imax=1e5, g=2)
transform <- rep("log",4)

## TESTING
obj(par_transform(estpars, transform), sdata, fixpars, transform)
traj_match(estpars, fixpars, transform, sdata, eval.onl=TRUE)$lik
pars <- c(estpars,fixpars)
mutate(sdata,
       E.F0=pars["F0"],
       E.Ft=calc_Ft(pars["a"],Z,pars["h"],pars["Imax"],pars["g"],pars["F0"],length,pars["fh"],t)) %>%
    with(.,
         dnorm(F0, mean=E.F0, sd=pars["Fobs"], log=TRUE) +
             dnorm(Ft, mean=E.Ft, sd=pars["Fobs"], log=TRUE)
         ) %>% sum

## Ranges for the estimated parameters
box <- cbind(lower=c(a=0.0001, h=1, Imax=1e2, g=1),
             upper=c(a=1, h=4, Imax=1e5, g=4))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=100000) %>%
    apply(., 1, as.list) %>%
        lapply(., unlist) -> guesses
## Pick a huge range of potential fh values, and perform the fitting
## estimating the other parameters
fh_vals <- seq(100, 10000, 100)
i <- 1
fixpars <- c(fh=fh_vals[i], F0=p[1], Fobs=p[2])
mclapply(guesses,traj_match,fixpars,transform,sdata,eval.only=TRUE,mc.cores=10) %>%
    lapply(., function(x) x$lik) %>% unlist -> guess_lik
guesses[order(guess_lik)[1:500]] -> refine
mclapply(refine,traj_match,fixpars,transform,sdata,method="Nelder-Mead",mc.cores=10) %>%
    lapply(., unlist) %>% unlist %>% matrix(., ncol=6, byrow=TRUE) %>% as.data.frame -> refine_lik
colnames(refine_lik) <- c(names(estpars),'lik','conv')
mutate(refine_lik, "log(a)"=log(a)) -> refine_lik
## these parameter values seem really strange. In particular, the parameters relating ingestion to spore load and host volume are all over the place. The parameter a, in particular, is all over the place, with estimates ranging from -400 to 400 on the log scale. Moreover, the parameter estimates don't make any sense: h > 3 would be crazy, as we are envisioning L^h to be something like body volume.
pairs(subset(refine_lik, lik < min(lik)+2)[,c(7,2:5)])

## However, the likelihood is better than when using the model that doesn't include the exp(-aZ/L^h) term.
bestp <- refine_lik[refine_lik$lik %>% which.min,] %>% unlist
with(sdata, exp(-bestp["a"]*Z/length^bestp["h"]))

## My guess is that a and h can trade off against one another, ensuring that exp(-a*Z/L^h) ends up at reasonable values.
## the fix, therefore, should be pretty simple: fix the value of h at a range of values (say 2 to 4) and estimate a only.
box <- cbind(lower=c(a=0.0001, Imax=1e2, g=1),
             upper=c(a=1, Imax=1e5, g=4))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=100000) %>%
    apply(., 1, as.list) %>%
        lapply(., unlist) -> guesses
fh_vals <- seq(100, 10000, 100)
h_vals <- seq(1,4,0.1)
est_params <- vector(mode='list', length=length(fh_vals))
for (i in 1:length(fh_vals)) {
    print(paste("fh =", fh_vals[i]))
    pars_vary_h <- vector(mode='list', length=length(h_vals))
    for (j in 1:length(h_vals)) {
        print(paste("h =", h_vals[j]))
        fixpars <- c(fh=fh_vals[i], F0=p[1], Fobs=p[2], h=h_vals[j])
        mclapply(guesses,traj_match,fixpars,transform,sdata,eval.only=TRUE,mc.cores=10) %>%
            lapply(., function(x) x$lik) %>% unlist -> guess_lik
        guesses[order(guess_lik)[1:500]] -> refine
        mclapply(refine,traj_match,fixpars,transform,sdata,method="Nelder-Mead",mc.cores=10) %>%
            lapply(., unlist) %>% unlist %>% matrix(., ncol=5, byrow=TRUE) %>% as.data.frame -> refine_lik
        colnames(refine_lik) <- c(names(guesses[[1]]),'lik','conv')
        pars_vary_h[[j]] <- refine_lik
    }
    est_params[[i]] <- pars_vary_h
    saveRDS(est_params, file="~/Dropbox/Infected_parameters_varying_h_and_fh.RDS")
}

## This seems to take care of the problem, except that both a and g are estimated to be very, very small numbers, indicating (again), that feeding is independent of


    mutate(refine_lik, "log(a)"=log(a)) -> refine_lik


est_params <- vector(mode='list', length=length(fh_vals))
for (i in 1:length(fh_vals)) {
    print(fixpars["fh"])
    t1 <- Sys.time()
    t2 <- Sys.time()
    print(paste("Guess lik calc", t2-t1))
    print(min(refine_lik$lik))
    t2 <- Sys.time()
    print(paste("Refine lik calc", t2-t1))
    est_params[[i]] <- refine_lik
    saveRDS(est_params, "Feeding_model_fitting_uninfecteds.RDS")
}
