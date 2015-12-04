source("Growth_reproduction_trajectory_fitting_functions_nondimensionalized_model.R")

## The parameters to be estimated (on the natural scale)
pars <- c(rho=0.1,
          K=0.5,
          km=0.01,
          ER=0.000225,
          EG=1.5,
          Lobs=0.1,
          W0=0.1)

## The fixed parameters
fh <- 1e-4/44.5e-9 ## about 2250 cells/ml
fixpars <- c(Imax=calc_Imax(fh), # fixed based on feeding data fitting and fh value
             fh=fh, # estimated from fitting growth/reproduction data
             g=calc_g(fh), # fixed based on feeding data fitting and fh value
             eps=44.5e-9, # fixed based on measured carbon content of algae
             V=30, # fixed based on experimental conditions
             xi=1.8e-3, # fixed based on Hall et al. 2009 length-weight regression
             q=3, # fixed based on Hall et al. 2009 length-weight regression
             v=100, # fixed based on the fact that it doesn't affect the fitting
             F0=1000000/30, # fixed based on experimental conditions
             L0=0.5, # fixed length at birth
             Wmat=0.002)

## The data
x <- read.csv("Cat_data/uninfected_growth_reproduction.csv")
data <- x[1:103,1:3]

## The feeding schedule
eventdat <- data.frame(var="F",
                       time=1:35,
                       value=unname(fixpars["F0"]),
                       method=rep(c(rep("add",4),"rep"),35/5))

## Testing the traj_match function
traj_match(pars=pars, obsdata=data, events=eventdat, eval.only=TRUE)

## Testing the obj function for confirmation
matrix(pars, nrow=1, byrow=TRUE, dimnames=list(NULL, names(pars))) %>%
        as.data.frame %>%
            with(.,
                 c(w_scalar=EG/K, ## reserve density scalar
                   l_scalar=unname(fixpars["v"])/km, ## length scalar
                   rho=rho, ## dimensionless growth efficiency
                   beta=EG/ER*(1-K)/K*unname(fixpars["xi"])*(unname(fixpars["v"])/km)^unname(fixpars["q"]), ## dimensionless birth rate
                   W0=W0, ## initial reserve density
                   Lobs=Lobs ## observation error for length
                   )
                 ) -> estpars
obj(estpars=par_transform(estpars,c("log","log","logit",rep("log",3))), data=data, fixpars=fixpars, events=eventdat)

## Note that the algorithm is passed a set of dimensional parameters, but it returns the dimensionless parameters. It is actually not possible to completely recover all of the dimensional parameters. In particular, km=v/l_scalar and rho=sigma*w_scalar. w_scalar is the ratio EG/K and beta/w_scalar is the ratio (1-K)/ER, but you cannot recover any of these parameter values on their own. This was actually indicated by the fitting as well, which suggested that ratios of K, EG, and ER might be well estimated, but not each parameter independently.

## initial guesses for the estimated parameters
box <- cbind(lower=c(rho=0, K=0, km=0.001, ER=1e-5, EG=0.5, Lobs=0.0001, W0=1e-3),
             upper=c(rho=1, K=1, km=10, ER=0.1, EG=5, Lobs=2, W0=1))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=500000) %>%
    apply(., 1, as.list) %>%
        lapply(., unlist) -> guesses

source("Growth_reproduction_trajectory_fitting_functions_nondimensionalized_model.R")
mclapply(guesses,
         traj_match,
         obsdata=data,
         events=eventdat,
         eval.only=TRUE,
         mc.cores=12) %>%
    lapply(., function(x) x$lik) %>%
        unlist -> guess_lik
guesses[order(guess_lik)[1:5000]] -> refine
mclapply(refine,
         traj_match,
         obsdata=data,
         events=eventdat,
         eval.only=FALSE,
         method="Nelder-Mead",
         mc.cores=12) -> refine_lik
refine_lik %>%
    lapply(., unlist) %>%
        unlist %>%
            matrix(., ncol=length(refine_lik[[1]]$params)+2, byrow=TRUE, dimnames=list(1:length(refine_lik), c(names(refine_lik[[1]]$params),"lik","conv"))) %>%
                as.data.frame -> refine_pars
refine_pars <- arrange(refine_pars, lik)

