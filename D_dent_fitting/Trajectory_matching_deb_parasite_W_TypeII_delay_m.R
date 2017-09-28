####################################################################
#################           RUN THIS          ######################
####################################################################
source("Growth_reproduction_trajectory_matching_infecteds.R")

## This model
model <- "deb_parasite_W_TypeII_delay_m"
dyn.load(paste0(model, ".so"))

x <- read.csv("Cat_data/all_data.csv") %>% as.data.frame
x2 <- subset(x, infected==1)
data <- x2[,c("day.of.assay","length..mm.","Spore.per.animal")]
colnames(data) <- c("age","length","spores")
## how to deal with 0s?
data$spores[is.na(data$spores)] <- 0

####################################################################
#################           TESTING           ######################
####################################################################
## parameters fixed at the uninfected data values
fixpars <- c(rho=0.152, K=1, v=10, F0=1e6/30, Lobs=0.102, P0=100)

## estimated parameters for this parasite model
estpars <- c(km=0.073, aP=1e-10, eP=1e10, hP=1e-2, m=0.01, tau=10, Pobs=500)

## unconstrained scale transformation
transform <- rep('log',7)

## ensure that the parameters are in the correct order for the model function
pars <- c(estpars, fixpars)
parorder <- c("rho","K","km","v","F0","Lobs","aP","eP","hP","m","tau","P0","Pobs")
pars[match(parorder, names(pars))] -> pars

## Generic set of food addition events
eventdat <- data.frame(var="F",
                       time=1:35,
                       value=unname(pars["F0"]),
                       method=rep(c(rep("add",4),"rep"),max(data$age)/5))

## Initial condition: assuming that E/W = rho/v and E+W=ER, then W
## + W*rho/v = ER, W(1+rho/v)=ER, W=ER/(1+rho/v)
y0 <- c(F=unname(pars["F0"]),
        E=0,
        W=unname(5.39e-5/(1+pars["rho"]/pars["v"])),
        Pi=unname(pars["P0"]),
        Pm=0)
y0["E"] <- unname(y0["W"]*pars["rho"]/pars["v"])

## Simulate the system
try(dede(y0,
        times=seq(0,35,0.1),
        func="derivs",
        parms=pars,
        dllname=model,
        initfunc="initmod",
        events=list(data=eventdat))) -> out
as.data.frame(out[out[,'time']%in%data$age,]) -> pred

## compute the observed weight as Wobs = W + E and compute the observed length prediction as Wobs=xi*Lobs^q; (Wobs/xi)^(1/q)=Lobs
xi <- 1.8e-3; q <- 3;
mutate(pred, Wobs=W+E, Lobs=(Wobs/xi)^(1/q)) -> pred

## compute the probability of observing the data, given the prediction
sapply(unique(data$age),
       function(d)
           c(dnorm(x=data$length[data$age==d],
                   mean=pred$Lobs[pred$time==d],
                   sd=pars["Lobs"],
                   log=TRUE) %>% sum,
             dnbinom(x=data$spores[data$age==d],
                     mu=pred$Pm[pred$time==d],
                     size=pars["Pobs"],
                     log=TRUE) %>% sum
             ) %>% sum
       ) %>% sum -> lik

tm_obj(par_transform(estpars, transform), data, fixpars, parorder, transform, model)

optimizer(estpars, fixpars, parorder, transform, data, model, eval.only=TRUE, method="Nelder-Mead")
tryout <- optimizer(estpars, fixpars, parorder, transform, data, model, eval.only=FALSE, method="Nelder-Mead")

####################################################################
#################           FITTING           ######################
####################################################################

## parameters fixed at the uninfected data values
fixpars <- c(rho=0.152, K=1, v=10, F0=1e6/30, Lobs=0.102, P0=100)

## parameter order expected by this model of parasite replication
parorder <- c("rho","K","km","v","F0","Lobs","aP","eP","hP","m","tau","P0", "Pobs")

## initial guesses at the estimated parameter values
box <- cbind(lower=c(km=0.073, aP=1e-12, eP=1e6, hP=1e-4, m=1e-4, tau=5, Pobs=10),
             upper=c(km=0.073, aP=1e-6, eP=1e12, hP=1e-1, m=1e-1, tau=15, Pobs=1000))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=600000) %>%
                apply(., 1, as.list) %>%
                    lapply(., unlist) -> guesses

## unconstrained scale transformation for estimated parameters
transform <- rep('log',7)

print("Estimating")
## Estimate likelihood at these guesses
mclapply(guesses,
         optimizer,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         model=model,
         eval.only=TRUE,
         type="trajectory_matching",
         mc.cores=12) %>%
             lapply(., function(x) x$lik) %>%
                 unlist -> guess_lik

print("Optimizing")
## optimize from the top 1000 parameter sets
guesses[order(guess_lik)[1:1000]] -> refine
mclapply(refine,
         optimizer,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         model=model,
         eval.only=FALSE,
         type='trajectory_matching',
         method='Nelder-Mead',
         mc.cores=15) -> refine_lik

saveRDS(refine_lik, file=paste0("TM_parameter_estimates_",model,".RDS"))

## parameters fixed at the uninfected data values
fixpars <- c(km=0.073, K=1, v=10, F0=1e6/30, Lobs=0.102, P0=100)

## parameter order expected by this model of parasite replication
parorder <- c("rho","K","km","v","F0","Lobs","aP","eP","m","tau","P0", "Pobs")

## initial guesses at the estimated parameter values
box <- cbind(lower=c(rho=0.152, aP=1e-12, eP=1e6, hP=1e-4, m=1e-4, tau=5, Pobs=10),
             upper=c(rho=0.152, aP=1e-6, eP=1e12, hP=1e-1, m=1e-1, tau=15, Pobs=1000))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=600000) %>%
                apply(., 1, as.list) %>%
                    lapply(., unlist) -> guesses

## unconstrained scale transformation for estimated parameters
transform <- c('logit',rep('log',6))

("Estimating")
## Estimate likelihood at these guesses
mclapply(guesses,
         optimizer,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         model=model,
         eval.only=TRUE,
         type="trajectory_matching",
         mc.cores=15) %>%
             lapply(., function(x) x$lik) %>%
                 unlist -> guess_lik

print("Optimizing")
## optimize from the top 1000 parameter sets
guesses[order(guess_lik)[1:1000]] -> refine
mclapply(refine,
         optimizer,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         model=model,
         eval.only=FALSE,
         type='trajectory_matching',
         method='Nelder-Mead',
         mc.cores=15) -> refine_lik

saveRDS(refine_lik, file=paste0("TM_parameter_estimates_",model,"_2.RDS"))

## parameters fixed at the uninfected data values
fixpars <- c(K=1, v=10, F0=1e6/30, Lobs=0.102, P0=100)

## parameter order expected by this model of parasite replication
parorder <- c("rho","K","km","v","F0","Lobs","aP","eP","m","tau","P0", "Pobs")

## initial guesses at the estimated parameter values
box <- cbind(lower=c(rho=0.152, km=0.073, aP=1e-12, eP=1e6, hP=1e-4, m=1e-4, tau=5, Pobs=10),
             upper=c(rho=0.152, km=0.073, aP=1e-6, eP=1e12, hP=1e-1, m=1e-1, tau=15, Pobs=1000))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=600000) %>%
                apply(., 1, as.list) %>%
                    lapply(., unlist) -> guesses

## unconstrained scale transformation for estimated parameters
transform <- c('logit',rep('log',7))

("Estimating")
## Estimate likelihood at these guesses
mclapply(guesses,
         optimizer,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         model=model,
         eval.only=TRUE,
         type="trajectory_matching",
         mc.cores=15) %>%
             lapply(., function(x) x$lik) %>%
                 unlist -> guess_lik

print("Optimizing")
## optimize from the top 1000 parameter sets
guesses[order(guess_lik)[1:1000]] -> refine
mclapply(refine,
         optimizer,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         model=model,
         eval.only=FALSE,
         type='trajectory_matching',
         method='Nelder-Mead',
         mc.cores=15) -> refine_lik

saveRDS(refine_lik, file=paste0("TM_parameter_estimates_",model,"_3.RDS"))
