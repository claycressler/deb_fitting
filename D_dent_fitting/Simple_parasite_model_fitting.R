source("Simple_parasite_model_fitting_functions.R")

system("rm Structured_parasite_model.so")
system("R CMD SHLIB Structured_parasite_model.c")

system("rm Unstructured_parasite_model.so")
system("R CMD SHLIB Unstructured_parasite_model.c")

## Simulate data
y0 <- c(E=100, C=4000, G=0, T=0) ## Initial conditions
## Baseline set of parameters - sample around these
pars <- c(theta=10, r=0.1, aC=0.0005, hC=20, aG=0.0001, hG=15, b=20000, m=10000, t_rep=15, obs_sd=100000)
datasets <- vector(mode='list', length=40)
i <- 1
set.seed(12340987)
while (is.null(datasets[[40]])) {
    this_pars <- rnorm(length(pars), mean=pars, sd=pars)
    while (any(this_pars < 0))
        this_pars <- rnorm(length(pars), mean=pars, sd=pars)
    names(this_pars) <- names(pars)
    this_pars["theta"] <- 10
    this_pars["t_rep"] <- round(this_pars["t_rep"])
    out <- ode(y0, times=seq(0, 50), func="derivs", parms=this_pars, dllname="Structured_parasite_model", initfunc="initmod", events=list(func="event", time=round(this_pars["t_rep"]))) ## simulated data
    if (all(out >= 0)) {
        data.frame(time=rep(seq(5,50,5), each=10),
                   spores=sapply(out[seq(5,50,5),"T"], function(m) rnorm(10, mean=m, sd=this_pars["obs_sd"])) %>% as.numeric) -> data
        data$spores[data$spores < 0] <- 0
        datasets[[i]] <- list(params=this_pars, data=data)
        i <- i+1
    }
}

str_parorder <- c("theta","r","aC","hC","aG","hG","b","m","t_rep","obs_sd")
unstr_parorder <- c("theta", "r", "aP", "hP", "b", "P0", "obs_sd")
unstr_fixpars <- str_fixpars <- c(theta=10)
str_transform <- rep("log",length(str_parorder)-1)
str_transform[which(str_parorder=="t_rep")-1] <- "logit"
unstr_transform <- rep("log",length(unstr_parorder)-1)

est_params <- vector(mode='list', length=40)
for (i in 1:40) {
    data <- datasets[[i]]$data

    ## STRUCTURED MODEL
    ## Generate a large number of different initial parameter guesses
    str_box <- cbind(lower=pars[-which(names(pars)=="theta")]/1000, upper=pars[-which(names(pars)=="theta")]*1000)
    str_box[rownames(str_box)=="t_rep",] <- c(0, 1)
    sobolDesign(lower=str_box[,"lower"], upper=str_box[,"upper"], nseq=500000) %>%
        apply(., 1, as.list) %>%
            lapply(., function(l) unlist(l)) -> guesses
    dyn.load("Structured_parasite_model.so")
    mclapply(guesses,
             traj_match,
             fixpars=str_fixpars,
             parorder=str_parorder,
             transform=str_transform,
             obsdata=data,
             eval.only=TRUE,
             method="subplex",
             mc.cores=10) %>%
        lapply(., function(x) x$lik) %>%
            unlist-> guess_lik
    ## Pick the best 1000 of these and use them as starting points for further refinement
    guesses[order(guess_lik)[1:1000]] -> refine
    mclapply(refine,
             traj_match,
             fixpars=str_fixpars,
             parorder=str_parorder,
             transform=str_transform,
             obsdata=data,
             eval.only=FALSE,
             method="Nelder-Mead",
             mc.cores=10) -> refine_lik
    refine_lik %>%
        lapply(., function(l) c(l$params, l$lik)) %>%
            unlist %>%
                matrix(., ncol=10, byrow=TRUE) %>%
                    as.data.frame -> refine_pars
    colnames(refine_pars) = c("r","aC","hC","aG","hG","b","m","t_rep","obs_sd", "lik")
    str_refine_pars <- arrange(refine_pars, lik)

    ## STRUCTURED MODEL
    ## Generate a large number of different initial parameter guesses
    unstr_box <- cbind(lower=c(r=1, aP=1e-4, hP=1, b=100000, P0=1000, obs_sd=1e5)/1000,
                       upper=c(r=1, aP=1e-4, hP=1, b=100000, P0=1000, obs_sd=1e5)*1000)
    sobolDesign(lower=unstr_box[,"lower"], upper=unstr_box[,"upper"], nseq=5000) %>%
        apply(., 1, as.list) %>%
            lapply(., function(l) unlist(l)) -> guesses
    ## compute the likelihood of each of these guesses
    dyn.load("Unstructured_parasite_model.so")
    mclapply(guesses,
             traj_match2,
             fixpars=unstr_fixpars,
             parorder=unstr_parorder,
             transform=unstr_transform,
             obsdata=data,
             eval.only=TRUE,
             method="Nelder-Mead",
             mc.cores=10) %>%
        lapply(., function(l) l$lik) %>%
            unlist -> guess_lik
    ## Pick the best 1000 for further refinement
    guesses[order(guess_lik)[1:1000]] -> refine
    mclapply(refine,
             traj_match2,
             fixpars=unstr_fixpars,
             parorder=unstr_parorder,
             transform=unstr_transform,
             obsdata=data,
             eval.only=FALSE,
             method="Nelder-Mead",
             mc.cores=10) -> refine_lik
    lapply(refine_lik, function(l) c(l$params, l$lik)) %>%
        unlist %>%
            matrix(., ncol=7, byrow=TRUE) %>%
                as.data.frame -> refine_pars
    colnames(refine_pars) = c( "r", "aP", "hP", "b", "P0", "obs_sd", "lik")
    arrange(refine_pars, lik) -> unstr_refine_pars

    est_params[[i]] <- list(str=str_refine_pars, unstr=unstr_refine_pars)
    saveRDS(est_params, file="Comparing_structred_unstructured_model_fits.RDS")
}
