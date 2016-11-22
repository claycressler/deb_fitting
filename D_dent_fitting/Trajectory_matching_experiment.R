source("Growth_reproduction_trajectory_matching.R")

if (!file.exists("Trajectory_matching_datasets_11-22.RDS")) {

    datasets <- vector(mode='list', length=25)
    set.seed(101)
    for (d in 1:25) {
        pars <- c(Imax=22500,
                  Fh=rnorm(1, mean=10000, sd=2000),
                  g=1.45,
                  rho=rnorm(1, mean=0.4, sd=0.08),
                  K=rnorm(1, mean=0.4, sd=0.1),
                  km=rnorm(1, mean=0.1, sd=0.03),
                  v=rnorm(1, mean=10, sd=1),
                  ER=rnorm(1, mean=0.001, sd=0.0005),
                  F0=1e6/30,
                  Lobs=0.1,
                  Robs=2,
                  Wmat=rnorm(1,mean=0.004,sd=0.0005))
        while(any(pars < 0))
            pars <- c(Imax=22500,
                      Fh=rnorm(1, mean=10000, sd=2000),
                      g=1.45,
                      rho=rnorm(1, mean=0.4, sd=0.08),
                      K=rnorm(1, mean=0.4, sd=0.1),
                      km=rnorm(1, mean=0.1, sd=0.03),
                      v=rnorm(1, mean=10, sd=1),
                      ER=rnorm(1, mean=0.001, sd=0.0005),
                      F0=1e6/30,
                      Lobs=0.1,
                      Robs=2,
                      Wmat=rnorm(1,mean=0.005,sd=0.001))
        pars["Imax"] <- calc_Imax(unname(pars["Fh"]))
        pars["g"] <- calc_g(unname(pars["Fh"]))
        ## days when observations will take place
        times <- c(5,10,12,15,18,25,30,35)
        ## number of reps per observation day
        inds <- 1:12
        ## Simulate 96 different trajectories to sample
        out <- vector(mode='list', length=length(times)*length(inds))
        for (i in 1:(length(times)*length(inds))) {
            ## feeding schedule - amount of food added each time is stochastic
            eventdat <- data.frame(var="F",
                                   time=1:35,
                                   value=rnorm(35,
                                       mean=unname(pars["F0"]),
                                       sd=5000),
                                   method=rep(c(rep("add",4),"rep"),7))
            y0 <- c(F=unname(pars["F0"]),
                    E=0,
                    W=unname(pars["ER"]/(1+pars["rho"]/pars["v"])),
                    R=0)
            y0["E"] <- unname(y0["W"]*pars["rho"]/pars["v"])
            ode(y0,
                times=seq(0,35,0.1),
                func="derivs",
                parms=pars,
                dllname="tm_deb",
                initfunc="initmod",
                events=list(data=eventdat)) -> out[[i]]
            ind <- which(out[[i]][,'W'] < pars['Wmat'])
            out[[i]][,'R'] <- out[[i]][,'R']-out[[i]][max(ind),'R']
            out[[i]][ind,'R'] <- 0
        }
        ## Simulate observations
        data = expand.grid(ind=inds, age=times)
        data$length <- rnorm(nrow(data),
                             mean=sapply(1:nrow(data),
                                 function(i)
                                     (sum(out[[i]][data[i,'age']*10+1,c('W','E')])/1.8e-3)^(1/3)
                                         ),
                             sd=pars['Lobs'])
        data$eggs <- sapply(1:nrow(data), function(i) out[[i]][data[i,'age']*10+1,'R'])
        data$eggs[data$eggs > 0] <- sapply(data$eggs[data$eggs > 0],
                                           function(e)
                                               rnorm(1, mean=e, sd=pars['Robs']))
        data$eggs[data$eggs < 0] <- 0
        datasets[[d]] <- list(data=data,
                              params=pars)
    }
    saveRDS(datasets, file="Trajectory_matching_datasets_11-22.RDS")
} else datasets <- readRDS("Trajectory_matching_datasets_11-22.RDS")

tm_ests <- vector(mode='list', length=25)
for (d in 1:25) {
    print(d)
    data <- datasets[[d]]$data

    ## Begin the funnel of optimization with simple trajectory matching,
    ## then use the particle filter to hone the parameter estimates.
    box <- cbind(lower=c(Fh=2000, rho=0, K=0, km=0.001, ER=0.00001, Lobs=0.001, Robs=0.01, Wmat=0.0001),
                 upper=c(Fh=20000, rho=1, K=1, km=1, ER=0.001, Lobs=2, Robs=10, Wmat=0.01))
    sobolDesign(lower=box[,'lower'],
                upper=box[,'upper'],
                nseq=50000) %>%
                    apply(., 1, as.list) %>%
                        lapply(., unlist) -> guesses

    fixpars <- c(Imax=22500, g=1.45, v=10, F0=1e6/30)
    transform <- c("log", rep("logit",2), rep("log",5))
    parorder <- c("Imax","Fh","g","rho","K","km","v","ER","F0","Lobs","Robs","Wmat")

    mclapply(guesses,
             optimizer,
             fixpars=fixpars,
             parorder=parorder,
             transform=transform,
             obsdata=data,
             eval.only=TRUE,
             type='trajectory_matching',
             mc.cores=15) %>%
                 lapply(., function(x) x$lik) %>%
                     unlist -> guess_lik
    guesses[order(guess_lik)[1:300]] -> refine
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
    tm_ests[[d]] <- refine_pars
    saveRDS(tm_ests, file="Trajectory_matching_estimates_11-22.RDS")
}

