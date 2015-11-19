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
lik

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

fixpars <- pars[c("Imax","g","eps","V","F0","xi","q","ER","v","Winit","Wmat")]
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
         mc.cores=4) %>%
        lapply(., function(x) x$lik) %>%
            unlist -> guess_lik
guesses[order(guess_lik)[1:5000]] -> refine
mclapply(refine,
         traj_match,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         events=eventdat,
         eval.only=FALSE,
         method="Nelder-Mead",
         mc.cores=4) -> refine_lik
refine_lik %>%
    lapply(., unlist) %>%
        unlist %>%
            matrix(., ncol=(nrow(box)+2), byrow=TRUE, dimnames=list(1:length(refine_lik), c(rownames(box),"lik","conv"))) %>%
                as.data.frame -> refine_pars
refine_pars <- arrange(refine_pars, lik)
saveRDS(refine_pars, file="~/Dropbox/Growth_reproduction_trajectory_fitting_Cat.RDS")
