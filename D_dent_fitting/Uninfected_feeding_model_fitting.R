source("Uninfected_feeding_model_fitting_functions.R")

data <- read.csv("Cat_data/feeding2.csv")
estpars <- c(Imax=1e5, g=2, F0=1e7, Fobs=1e4)
fixpars <- c(fh=2.5e6)
transform <- rep("log",4)
## TESTING
obj(par_transform(estpars, transform), data, fixpars, transform)
traj_match(estpars, fixpars, transform, data, eval.onl=TRUE)
pars <- c(estpars,fixpars)
mutate(data,
       E.F0=estpars["F0"],
       E.Ft=calc_Ft(pars["Imax"],pars["g"],pars["F0"],length,pars["fh"],t)) %>%
    with(.,
         dnorm(F0, mean=E.F0, sd=pars["Fobs"], log=TRUE) +
             dnorm(Ft, mean=E.Ft, sd=pars["Fobs"], log=TRUE)
         ) %>% sum

## Ranges for the estimated parameters
box <- cbind(lower=c(Imax=1e2, g=1, F0=5000, Fobs=1e2),
             upper=c(Imax=1e5, g=4, F0=15000, Fobs=1e3))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=100000) %>%
    apply(., 1, as.list) %>%
        lapply(., unlist) -> guesses

run <- FALSE
if (run) {
## Fitting the feeding parameters of uninfected animals
data <- subset(data, infected==0)
## Pick a huge range of potential fh values, and perform the fitting
## estimating the other parameters
fh_vals <- seq(100, 10000, 100)
est_params <- vector(mode='list', length=length(fh_vals))
for (i in 1:length(fh_vals)) {
    fixpars <- c(fh=fh_vals[i])
    print(fixpars["fh"])
    t1 <- Sys.time()
    mclapply(guesses,
             traj_match,
             fixpars,
             transform,
             data,
             eval.only=TRUE,
             mc.cores=5) %>%
        lapply(., function(x) x$lik) %>%
            unlist -> guess_lik
    t2 <- Sys.time()
    print(paste("Guess lik calc", t2-t1))
    guesses[order(guess_lik)[1:500]] -> refine
    t1 <- Sys.time()
    mclapply(refine,
             traj_match,
             fixpars,
             transform,
             data,
             eval.only=FALSE,
             method="Nelder-Mead",
             mc.cores=5) %>%
        lapply(., unlist) %>%
            unlist %>%
                matrix(., ncol=6, byrow=TRUE) %>%
                    as.data.frame -> refine_lik
    colnames(refine_lik) <- c(names(estpars),'lik','conv')
    print(min(refine_lik$lik))
    t2 <- Sys.time()
    print(paste("Refine lik calc", t2-t1))
    est_params[[i]] <- refine_lik
    saveRDS(est_params, "Feeding_model_fitting_uninfecteds.RDS")
}

## Fitting the feeding parameters of infected animals
data <- subset(data, infected==1)
## Pick a huge range of potential fh values, and perform the fitting
## estimating the other parameters
fh_vals <- seq(100, 10000, 100)
est_params <- vector(mode='list', length=length(fh_vals))
for (i in 1:length(fh_vals)) {
    fixpars <- c(fh=fh_vals[i])
    print(fixpars["fh"])
    t1 <- Sys.time()
    mclapply(guesses,
             traj_match,
             fixpars,
             transform,
             data,
             eval.only=TRUE,
             mc.cores=10) %>%
        lapply(., function(x) x$lik) %>%
            unlist -> guess_lik
    t2 <- Sys.time()
    print(paste("Guess lik calc", t2-t1))
    guesses[order(guess_lik)[1:500]] -> refine
    t1 <- Sys.time()
    mclapply(refine,
             traj_match,
             fixpars,
             transform,
             data,
             eval.only=FALSE,
             method="Nelder-Mead",
             mc.cores=10) %>%
        lapply(., unlist) %>%
            unlist %>%
                matrix(., ncol=6, byrow=TRUE) %>%
                    as.data.frame -> refine_lik
    colnames(refine_lik) <- c(names(estpars),'lik','conv')
    print(min(refine_lik$lik))
    t2 <- Sys.time()
    print(paste("Refine lik calc", t2-t1))
    est_params[[i]] <- refine_lik
    saveRDS(est_params, "Feeding_model_fitting_infecteds_all_pars.RDS")
}
}

## Fitting the feeding parameters of infected animals, assuming only
## Imax differs between infected and uninfected animals
data <- read.csv("Cat_data/feeding2.csv")
data <- subset(data, infected==1)
fh_vals <- seq(100, 10000, 100)
x <- readRDS("Feeding_model_fitting_uninfecteds.RDS")
sapply(1:length(fh_vals),
       function(i) {
           fixpars <- c(fh=fh_vals[i],
                        subset(x[[i]], lik==min(lik))[c("g","F0","Fobs")]) %>% unlist
           traj_match(c(Imax=7500),
                      fixpars,
                      transform=c('log'),
                      obsdata=data) %>%
               unlist
       }) -> inf_Imax_est
saveRDS(inf_Imax_est, "Feeding_model_fitting_infecteds_Imax_only.RDS")
