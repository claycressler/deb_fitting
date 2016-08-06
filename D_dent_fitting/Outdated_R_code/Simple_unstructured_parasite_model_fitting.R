source("Simple_parasite_model_fitting_functions.R")

system("rm Unstructured_parasite_model.so")
system("R CMD SHLIB Unstructured_parasite_model.c")
dyn.load("Unstructured_parasite_model.so")
fixpars <- c(theta=10)
estpars <- c(r=1, aP=1e-4, hP=1, b=100000, P0=1000, obs_sd=1e5)
parorder <- c("theta", "r", "aP", "hP", "b", "P0", "obs_sd")
transform <- rep("log",6)
## testing
all(estpars-par_untransform(par_transform(estpars, transform), transform)<1e-10)
unstr_obj(par_transform(estpars,transform), data=data, fixpars=fixpars, parorder=parorder, transform=transform)==traj_match2(estpars, fixpars, parorder, transform, data, eval.only=TRUE)$lik

## Generate a ton of parameter estimates
library(pomp)
box <- cbind(lower=estpars/1000, upper=estpars*1000)
sobolDesign(lower=box[,"lower"], upper=box[,"upper"], nseq=500000) %>%
    apply(., 1, as.list) %>%
        lapply(., function(l) unlist(l)) -> guesses
## compute the likelihood of each of these guesses
mclapply(guesses,
         traj_match2,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         eval.only=TRUE,
         method="Nelder-Mead",
         mc.cores=10) -> guess_lik
guess_lik %>%
    lapply(., function(l) l$lik) %>%
        unlist -> guess_lik

## Pick the best 1000 for further refinement
guesses[order(guess_lik)[1:1000]] -> refine
mclapply(refine,
         traj_match2,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         eval.only=FALSE,
         method="Nelder-Mead",
         mc.cores=10) -> refine_lik

lapply(refine_lik, function(l) c(l$params, l$lik)) %>%
    unlist %>%
        matrix(., ncol=7, byrow=TRUE) %>%
            as.data.frame -> refine_pars
colnames(refine_pars) = c(names(estpars), "lik")
arrange(refine_pars, lik) -> refine_pars
saveRDS(refine_pars, file="Simulated_data_fits_Untructured_parasite_model.RDS")
