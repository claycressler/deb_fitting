source("Simple_parasite_model_fitting_functions.R")

system("rm Structured_parasite_model.so")
system("R CMD SHLIB Structured_parasite_model.c")
dyn.load("Structured_parasite_model.so")
## Simulate data
y0 <- c(E=100, C=4000, G=0, T=0) ## Initial conditions
str_params <- c(theta=10, r=0.1, aC=0.00005, hC=1, aG=0.00001, hG=0.5, b=200000, m=10000, t_rep=15, obs_sd=100000) ## True parameter values
out <- ode(y0, times=seq(0, 50), func="derivs", parms=str_params, dllname="Structured_parasite_model", initfunc="initmod", events=list(func="event", time=str_params["t_rep"])) ## simulated data
set.seed(100011011)
data.frame(time=rep(seq(5,50,5), each=10),
           spores=sapply(out[seq(5,50,5),"T"], function(m) rnorm(10, mean=m, sd=str_params["obs_sd"])) %>% as.numeric) -> data ## simulated observations
data$spores[data$spores<0] <- 0

## I will estimate all of the parameters *except* theta, which is the
## rate that energy becomes available to the system, and is thus
## analogous to ingestion in the full DEB model. The parameters of the
## ingestion model are estimated separately.
parorder <- c("theta","r","aC","hC","aG","hG","b","m","t_rep","obs_sd")
estpars <- str_params[-(((str_params %>% names)=="theta") %>% which)]
estpars["t_rep"] <- estpars["t_rep"]/max(data$time)
fixpars <- str_params["theta"]
transform <- rep("log",length(estpars))
transform[which(names(estpars)=="t_rep")] <- "logit"
## Testing
all(estpars-par_untransform(par_transform(estpars, transform), transform) < 1e-10)
str_obj(par_transform(estpars, transform), data, fixpars, parorder, transform)==
    traj_match(estpars=estpars, fixpars=fixpars, parorder=parorder, transform=transform, obsdata=data, eval.only=TRUE, method="subplex")$lik

## Generate a large number of different initial guesses distributed around the true parameter values
## Generate a ton of parameter estimates
box <- cbind(lower=estpars/100, upper=estpars*100)
sobolDesign(lower=box[,"lower"], upper=box[,"upper"], nseq=500000) %>%
    apply(., 1, as.list) %>%
        lapply(., function(l) unlist(l)) -> guesses
mclapply(guesses,
         traj_match,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         eval.only=TRUE,
         method="subplex",
         mc.cores=10) %>%
    lapply(., function(x) x$lik) %>%
        unlist-> guess_lik
## Most of these guesses were total shit
## Pick the best 1000 of these and use them as starting points for further refinement
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
    lapply(., function(l) c(l$params, l$lik)) %>%
        unlist %>%
            matrix(., ncol=10, byrow=TRUE) %>%
                as.data.frame -> refine_pars
colnames(refine_pars) = c(names(estpars), "lik")
refine_pars <- arrange(refine_pars, lik)
saveRDS(refine_pars, file="Simulated_data_fits_Structured_parasite_model.RDS")
