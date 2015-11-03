source("Infected_feeding_model_fitting_functions.R")

## FIT ONLY INDIVIDUALS WITH NON-NA SPORE COUNTS
data <- read.csv("Cat_data/feeding_infection.csv")
sdata <- subset(data, infected==1 & !is.na(Z))

## extract some parameters from previous fits to remove those from the estimation procedure
readRDS("Feeding_model_fitting_uninfecteds.RDS") %>%
    lapply(., function(l) subset(l, lik==min(lik))[c("g","F0","Fobs")]) %>%
        unlist %>% matrix(., ncol=3, byrow=T) %>% apply(., 2, mean) -> p
fixpars <- c(fh=1000, F0=p[2], Fobs=p[3], g=p[1])
estpars <- c(a=1, h=2, Imax=1e5)
transform <- rep("log",3)

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
box <- cbind(lower=c(a=0.0001, h=1, Imax=1e2),
             upper=c(a=1, h=4, Imax=1e5))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=100000) %>%
    apply(., 1, as.list) %>%
        lapply(., unlist) -> guesses
## Pick a huge range of potential fh values, and perform the fitting
## estimating the other parameters
fh_vals <- seq(100, 10000, 100)
i <- 1
fixpars <- c(fh=fh_vals[i], g=p[1], F0=p[2], Fobs=p[3])
mclapply(guesses,traj_match,fixpars,transform,sdata,eval.only=TRUE,mc.cores=10) %>%
    lapply(., function(x) x$lik) %>% unlist -> guess_lik
guesses[order(guess_lik)[1:500]] -> refine
mclapply(refine,traj_match,fixpars,transform,sdata,method="Nelder-Mead",mc.cores=10) %>%
    lapply(., unlist) %>% unlist %>% matrix(., ncol=5, byrow=TRUE) %>% as.data.frame -> refine_lik
colnames(refine_lik) <- c(names(estpars),'lik','conv')
mutate(refine_lik, "log(a)"=log(a)) -> refine_lik
pairs(subset(refine_lik, lik < min(lik)+1)[,c(1:4)])

est_params <- vector(mode='list', length=length(fh_vals))
for (i in 1:length(fh_vals)) {
    print(paste("fh =", fh_vals[i]))
    t1 <- Sys.time()
    ## extract the best fitting values of g, F0, and Fobs for this value of fh
    readRDS("Feeding_model_fitting_uninfecteds.RDS")[[i]] %>%
        subset(., lik < min(lik)+0.1) %>% ## extract all parameter sets near max lik
            apply(., 2, mean) -> p ## take the mean
    fixpars <- c(fh=fh_vals[i], p["g"], p["F0"], p["Fobs"])
    mclapply(guesses,traj_match,fixpars,transform,sdata,eval.only=TRUE,mc.cores=10) %>%
        lapply(., function(x) x$lik) %>% unlist -> guess_lik
    guesses[order(guess_lik)[1:500]] -> refine
    mclapply(refine,traj_match,fixpars,transform,sdata,method="Nelder-Mead",mc.cores=10) %>%
        lapply(., unlist) %>% unlist %>% matrix(., ncol=5, byrow=TRUE) %>% as.data.frame -> refine_lik
    colnames(refine_lik) <- c(names(guesses[[1]]),'lik','conv')
    est_params[[i]] <- refine_lik
    saveRDS(est_params, file="Feeding_model_fitting_infecteds_spore_dependent.RDS")
    t2 <- Sys.time()
    print(t2-t1)
}


## FIT ALL INDIVIDUALS - SET NA SPORE COUNTS TO 0
data <- read.csv("Cat_data/feeding_infection.csv")
sdata <- subset(data, infected==1)
sdata$Z[is.na(sdata$Z)] = 0

## Ranges for the estimated parameters
box <- cbind(lower=c(a=0.0001, h=1, Imax=1e2),
             upper=c(a=1, h=4, Imax=1e5))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=100000) %>%
    apply(., 1, as.list) %>%
        lapply(., unlist) -> guesses
transform <- rep("log",3)
## Pick a huge range of potential fh values, and perform the fitting
## estimating the other parameters
fh_vals <- seq(100, 10000, 100)
est_params <- vector(mode='list', length=length(fh_vals))
for (i in 1:length(fh_vals)) {
    print(paste("fh =", fh_vals[i]))
    t1 <- Sys.time()
    ## extract the best fitting values of g, F0, and Fobs for this value of fh
    readRDS("Feeding_model_fitting_uninfecteds.RDS")[[i]] %>%
        subset(., lik < min(lik)+0.1) %>% ## extract all parameter sets near max lik
            apply(., 2, mean) -> p ## take the mean
    fixpars <- c(fh=fh_vals[i], p["g"], p["F0"], p["Fobs"])
    mclapply(guesses,traj_match,fixpars,transform,sdata,eval.only=TRUE,mc.cores=10) %>%
        lapply(., function(x) x$lik) %>% unlist -> guess_lik
    guesses[order(guess_lik)[1:500]] -> refine
    mclapply(refine,traj_match,fixpars,transform,sdata,method="Nelder-Mead",mc.cores=10) %>%
        lapply(., unlist) %>% unlist %>% matrix(., ncol=5, byrow=TRUE) %>% as.data.frame -> refine_lik
    colnames(refine_lik) <- c(names(guesses[[1]]),'lik','conv')
    est_params[[i]] <- refine_lik
    saveRDS(est_params, file="Feeding_model_fitting_infecteds_spore_dependent_2.RDS")
    t2 <- Sys.time()
    print(t2-t1)
}
