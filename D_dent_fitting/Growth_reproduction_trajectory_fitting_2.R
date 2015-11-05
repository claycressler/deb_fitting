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

