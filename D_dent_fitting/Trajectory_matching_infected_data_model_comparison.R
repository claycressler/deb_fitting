## Energy from reserves
## Fitting aP only
source("Growth_reproduction_trajectory_matching_infecteds_model_1v2.R")

x <- read.csv("Cat_data/all_data.csv") %>% as.data.frame
x2 <- subset(x, infected==1)
data <- x2[,c("day.of.assay","length..mm.","Spore.per.animal")]
colnames(data) <- c("age","length","spores")
## how to deal with 0s?
data$spores[is.na(data$spores)] <- 0

## what if you fix everything except aP?
fixpars <- c(rho=0.152, K=1, km=0.073, v=10, F0=1e6/30, Lobs=0.0102, eP=1/(0.48*33e-9), Pobs=1000)
## unconstrained scale transformation
transform <- 'log'
parorder <- c("rho","K","km","v","F0","Lobs","aP","eP","Pobs")

## Estimate the parameters
guesses = as.list(seq(1e-10, 1e-2, length=10000))
for (i in 1:10000) names(guesses[[i]]) <- "aP"

## Estimate likelihood
mclapply(guesses,
         optimizer,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         eval.only=TRUE,
         type="trajectory_matching",
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
         method='subplex',
         mc.cores=15) -> refine_lik
saveRDS(refine_lik, file="TM_estimates_parasite_model_1_aP_only.RDS")

## min -log-likelihood: 649,067.5 > 430,409.9

pars <- c(refine_lik[[1]]$params, fixpars)
pars[match(parorder, names(pars))] -> pars

y0 <- c(F=unname(pars["F0"]),
        E=0,
        W=unname(5.39e-5/(1+pars["rho"]/pars["v"])),
        P=10)
y0["E"] <- unname(y0["W"]*pars["rho"]/pars["v"])

try(ode(y0,
        times=seq(0,35,0.1),
        func="derivs",
        parms=pars,
        dllname="tm_deb_parasite_1v2",
        initfunc="initmod",
        events=list(data=eventdat))) %>% as.data.frame %>%
            mutate(., Wobs=W+E, Lobs=(Wobs/1.8e-3)^(1/3)) -> out

## What if you allow parasite manipulation of somatic maintenance rate?
fixpars <- c(rho=0.152, K=1, v=10, F0=1e6/30, Lobs=0.0102, eP=1/(0.48*33e-9), Pobs=1000)
## unconstrained scale transformation
transform <- rep('log',2)
parorder <- c("rho","K","km","v","F0","Lobs","aP","eP","Pobs")

box <- cbind(lower=c(km=0.001, aP=1e-9),
             upper=c(km=1, aP=1e-4))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=100000) %>%
                apply(., 1, as.list) %>%
                    lapply(., unlist) -> guesses

## Estimate likelihood
mclapply(guesses,
         optimizer,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         eval.only=TRUE,
         type="trajectory_matching",
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
saveRDS(refine_lik, file="TM_estimates_parasite_model_1_aP_km.RDS")

####################################################################3
## Energy from growth allocation
## fitting phi only
source("Growth_reproduction_trajectory_matching_infecteds_model_2.R")

x <- read.csv("Cat_data/all_data.csv") %>% as.data.frame
x2 <- subset(x, infected==1)
data <- x2[,c("day.of.assay","length..mm.","Spore.per.animal")]
colnames(data) <- c("age","length","spores")
## how to deal with 0s?
data$spores[is.na(data$spores)] <- 0

## what if you fix everything except aP?
fixpars <- c(rho=0.152, K=1, km=0.073, v=10, F0=1e6/30, Lobs=0.0102, eP=1/(0.48*33e-9), Pobs=1000)
## unconstrained scale transformation
transform <- 'logit'
parorder <- c("rho","K","km","v","F0","Lobs","phi","eP","Pobs")

## Estimate the parameters
guesses = as.list(seq(0, 0.99999, length=10000))
for (i in 1:10000) names(guesses[[i]]) <- "phi"

## Estimate likelihood
mclapply(guesses,
         optimizer,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         eval.only=TRUE,
         type="trajectory_matching",
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
         method='subplex',
         mc.cores=15) -> refine_lik
saveRDS(refine_lik, file="TM_estimates_parasite_model_2_phi_only.RDS")


####################################################################3
## Energy from soma
## Fitting aP only
source("Growth_reproduction_trajectory_matching_infecteds_model_3.R")

x <- read.csv("Cat_data/all_data.csv") %>% as.data.frame
x2 <- subset(x, infected==1)
data <- x2[,c("day.of.assay","length..mm.","Spore.per.animal")]
colnames(data) <- c("age","length","spores")
## how to deal with 0s?
data$spores[is.na(data$spores)] <- 0

## what if you fix everything except aP?
fixpars <- c(rho=0.152, K=1, km=0.073, v=10, F0=1e6/30, Lobs=0.0102, eP=1/(0.48*33e-9), Pobs=5000)
## unconstrained scale transformation
transform <- rep('log',1)
parorder <- c("rho","K","km","v","F0","Lobs","aP","eP","Pobs")

## Estimate the parameters
guesses = as.list(seq(1e-10, 1e-2, length=10000))
for (i in 1:10000) names(guesses[[i]]) <- "aP"

## Estimate likelihood
mclapply(guesses,
         optimizer,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         eval.only=TRUE,
         type="trajectory_matching",
         mc.cores=15) %>%
             lapply(., function(x) x$lik) %>%
                 unlist -> guess_lik
guesses[order(guess_lik)[1:30]] -> refine
mclapply(refine,
         optimizer,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         eval.only=FALSE,
         type='trajectory_matching',
         method='subplex',
         mc.cores=15) -> refine_lik
saveRDS(refine_lik, file="TM_estimates_parasite_model_3_aP_only.RDS")

## min -log-likelihood: 430409.9

pars <- c(refine_lik[[1]]$params, fixpars)
pars[match(parorder, names(pars))] -> pars

y0 <- c(F=unname(pars["F0"]),
        E=0,
        W=unname(5.39e-5/(1+pars["rho"]/pars["v"])),
        P=10)
y0["E"] <- unname(y0["W"]*pars["rho"]/pars["v"])

try(ode(y0,
        times=seq(0,35,0.1),
        func="derivs",
        parms=pars,
        dllname="tm_deb_parasite_3",
        initfunc="initmod",
        events=list(data=eventdat))) %>% as.data.frame %>%
            mutate(., Wobs=W+E, Lobs=(Wobs/1.8e-3)^(1/3)) -> out

## What if you allow parasite manipulation of somatic maintenance rate?
fixpars <- c(rho=0.152, K=1, v=10, F0=1e6/30, Lobs=0.0102, eP=1/(0.48*33e-9), Pobs=5000)
## unconstrained scale transformation
transform <- rep('log',2)
parorder <- c("rho","K","km","v","F0","Lobs","aP","eP","Pobs")

box <- cbind(lower=c(km=0.001, aP=1e-9),
             upper=c(km=1, aP=1e-4))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=50000) %>%
                apply(., 1, as.list) %>%
                    lapply(., unlist) -> guesses

## Estimate likelihood
mclapply(guesses,
         optimizer,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         eval.only=TRUE,
         type="trajectory_matching",
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
saveRDS(refine_lik, file="TM_estimates_parasite_model_3_aP_km.RDS")
## min -log-likelihood: 420104.6 < 430409.9 - this version is better!

pars <- c(refine_lik[[1]]$params, fixpars)
pars[match(parorder, names(pars))] -> pars

y0 <- c(F=unname(pars["F0"]),
        E=0,
        W=unname(5.39e-5/(1+pars["rho"]/pars["v"])),
        P=10)
y0["E"] <- unname(y0["W"]*pars["rho"]/pars["v"])

try(ode(y0,
        times=seq(0,35,0.1),
        func="derivs",
        parms=pars,
        dllname="tm_deb_parasite_3",
        initfunc="initmod",
        events=list(data=eventdat))) %>% as.data.frame %>%
            mutate(., Wobs=W+E, Lobs=(Wobs/1.8e-3)^(1/3)) -> out


####################################################################3
## Energy from soma
## saturating functional response
## Fitting aP only
source("Growth_reproduction_trajectory_matching_infecteds_model_3v2.R")

x <- read.csv("Cat_data/all_data.csv") %>% as.data.frame
x2 <- subset(x, infected==1)
data <- x2[,c("day.of.assay","length..mm.","Spore.per.animal")]
colnames(data) <- c("age","length","spores")
## how to deal with 0s?
data$spores[is.na(data$spores)] <- 0

## what if you fix everything except aP?
fixpars <- c(rho=0.152, K=1, km=0.073, v=10, F0=1e6/30, Lobs=0.0102, eP=1/(0.48*33e-9), Pobs=1000)
## unconstrained scale transformation
transform <- rep('log',2)
parorder <- c("rho","K","km","v","F0","Lobs","aP","eP","Pobs","hP")

box <- cbind(lower=c(aP=1e-9, hP=1e-5),
             upper=c(aP=1e-4, hP=1e-2))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=50000) %>%
                apply(., 1, as.list) %>%
                    lapply(., unlist) -> guesses

## Estimate likelihood
mclapply(guesses,
         optimizer,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         eval.only=TRUE,
         type="trajectory_matching",
         mc.cores=15) %>%
             lapply(., function(x) x$lik) %>%
                 unlist -> guess_lik
guesses[order(guess_lik)[1:30]] -> refine
mclapply(refine,
         optimizer,
         fixpars=fixpars,
         parorder=parorder,
         transform=transform,
         obsdata=data,
         eval.only=FALSE,
         type='trajectory_matching',
         method='subplex',
         mc.cores=15) -> refine_lik
saveRDS(refine_lik, file="TM_estimates_parasite_model_3v2_aP_only.RDS")

## min -log-likelihood: 430409.9

pars <- c(refine_lik[[30]]$params, fixpars)
pars[match(parorder, names(pars))] -> pars

y0 <- c(F=unname(pars["F0"]),
        E=0,
        W=unname(5.39e-5/(1+pars["rho"]/pars["v"])),
        P=10)
y0["E"] <- unname(y0["W"]*pars["rho"]/pars["v"])

try(ode(y0,
        times=seq(0,35,0.1),
        func="derivs",
        parms=pars,
        dllname="tm_deb_parasite_3v2",
        initfunc="initmod",
        events=list(data=eventdat))) %>% as.data.frame %>%
            mutate(., Wobs=W+E, Lobs=(Wobs/1.8e-3)^(1/3)) -> out
out[out$time%in%data$age,] -> pred


##############################################

## Make some plots baby!

## Real data
x <- read.csv("Cat_data/all_data.csv") %>% as.data.frame
x2 <- subset(x, infected==1)
data <- x2[,c("day.of.assay","length..mm.","Spore.per.animal")]
colnames(data) <- c("age","length","spores")
## how to deal with 0s?
data$spores[is.na(data$spores)] <- 0

## Plot of feeding data
library(ggplot2)
ggplot(x, aes(x=day.of.assay, y=Clearance.rate, group=as.factor(infected), col=as.factor(infected))) +
    geom_point() +
        stat_smooth() +
            theme_bw() +
                theme(legend.position="none") +
                    xlab("Age") + ylab("Feeding rate") +
                        scale_colour_manual(values=c("black","red"))

ggplot(data, aes(x=age, y=spores)) +
    geom_point() +
        theme_bw() +
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=16)) +
                      xlab("Age") + ylab("Parasite burden")

## Feeding model
data <- read.csv("Cat_data/feeding_infection.csv")
## focus on uninfected animals first
uninf_data <- subset(data, infected==0)
## calculate the expected F0 and Ft
source("Uninfected_feeding_model_fitting_functions.R")
p0 <- readRDS("RDS_files/Feeding_model_fitting_uninfecteds.RDS")[[50]] %>% subset(., lik==min(lik))
pars <- c(p0[1:4], fh=5000) %>% unlist
## calculate the clearance rate as Cat did
mutate(uninf_data,
       E.F0=pars["F0"],
       E.Ft=calc_Ft(pars["Imax"],pars["g"],pars["F0"],length,pars["fh"],t),
       clearance=log(F0/Ft)*0.01/0.125,
       E.clearance=log(E.F0/E.Ft)*0.01/0.125) -> uninf_data
mutate(uninf_data, diff=clearance-E.clearance) -> uninf_data

with(uninf_data, plot(day, clearance, pch=21, bg=1))
with(inf_data, points(day, clearance, pch=21, bg=2))

for (i in 1:8)
    lines(x=rep(ddply(uninf_data, .(day), summarize, min.E=min(E.clearance), max.E=max(E.clearance))[i,1],2),
          y=ddply(uninf_data, .(day), summarize, min.E=min(E.clearance), max.E=max(E.clearance))[i,2:3],
          lwd=2)


## change spore load to zero for all infected animals, even if they weren't counted
inf_data <- subset(data, infected==1)
inf_data$Z[is.na(inf_data$Z)] <- 0
## calculate the expected F0 and Ft
source("Infected_feeding_model_fitting_functions.R")
p1 <- readRDS("RDS_files/Feeding_model_fitting_infecteds_spore_dependent_2.RDS")[[50]] %>% subset(., lik==min(lik))
pars1 <- c(pars["g"],pars["F0"],pars["Fobs"],pars["fh"],p1[1:3]) %>% unlist
## calculate the clearance rate as Cat did
mutate(inf_data,
       E.F0=pars1["F0"],
       E.Ft=calc_Ft(pars1["a"],Z,pars1["h"],pars1["Imax"],pars1["g"],pars1["F0"],length,pars["fh"],t),
       clearance=log(F0/Ft)*0.01/0.125,
       E.clearance=log(E.F0/E.Ft)*0.01/0.125) -> inf_data




rbind(gather(uninf_data, clear, rate, (ncol(uninf_data)-1):ncol(uninf_data)),
      gather(inf_data, clear, rate, (ncol(inf_data)-1):ncol(inf_data))) -> dat
mutate(dat, type=paste0(clear,infected)) -> dat

ggplot(dat, aes(x=day, y=rate, shape=clear, colour=factor(infected))) +
    geom_point() +
        theme_bw() +
            scale_colour_manual(values=c("black","red"), name="", labels=c("Uninfected","Infected")) +
                scale_shape_discrete(solid=F, name="", labels=c("Observed","Predicted")) +
                    xlab("Age") + ylab("Clearance rate")



## Best model, based on likelihoods
x <- read.csv("Cat_data/all_data.csv") %>% as.data.frame
x2 <- subset(x, infected==1)
data <- x2[,c("day.of.assay","length..mm.","Spore.per.animal")]
colnames(data) <- c("age","length","spores")
## how to deal with 0s?
data$spores[is.na(data$spores)] <- 0

x1 <- readRDS("TM_estimates_parasite_model_3_aP_only.RDS")
x2 <- readRDS("TM_estimates_parasite_model_3_aP_km.RDS")

fixpars <- c(rho=0.152, K=1, km=0.073, v=10, F0=1e6/30, Lobs=0.0102, eP=1/(0.48*33e-9), Pobs=5000)
fixpars1 <- c(rho=0.152, K=1, v=10, F0=1e6/30, Lobs=0.0102, eP=1/(0.48*33e-9), Pobs=5000)
parorder <- c("rho","K","km","v","F0","Lobs","aP","eP","Pobs")

pars <- c(x1[[lapply(x1, function(l) l$lik) %>% unlist %>% which.min]]$params, fixpars)
pars[match(parorder, names(pars))] -> pars
pars1 <- c(x2[[lapply(x2, function(l) l$lik) %>% unlist %>% which.min]]$params, fixpars1)
pars1[match(parorder, names(pars1))] -> pars1

y0 <- c(F=unname(pars["F0"]),
        E=0,
        W=unname(5.39e-5/(1+pars["rho"]/pars["v"])),
        P=10)
y0["E"] <- unname(y0["W"]*pars["rho"]/pars["v"])

try(ode(y0,
        times=seq(0,35,0.1),
        func="derivs",
        parms=pars,
        dllname="tm_deb_parasite_3",
        initfunc="initmod",
        events=list(data=eventdat))) %>% as.data.frame %>%
            mutate(., Wobs=W+E, Lobs=(Wobs/1.8e-3)^(1/3)) -> out

try(ode(y0,
        times=seq(0,35,0.1),
        func="derivs",
        parms=pars1,
        dllname="tm_deb_parasite_3",
        initfunc="initmod",
        events=list(data=eventdat))) %>% as.data.frame %>%
            mutate(., Wobs=W+E, Lobs=(Wobs/1.8e-3)^(1/3)) -> out1


## comparison of growth and spore production to actual data
bootL <- sapply(out[unique(data$age)*10+1,"Lobs"], function(m) rnorm(1000, mean=m, sd=pars["Lobs"]))
bootP <- sapply(out[unique(data$age)*10+1,"P"], function(m) rgamma(1000, shape=m/pars["Pobs"], scale=pars["Pobs"]))

bootL1 <- sapply(out1[unique(data$age)*10+1,"Lobs"], function(m) rnorm(1000, mean=m, sd=pars1["Lobs"]))
bootP1 <- sapply(out1[unique(data$age)*10+1,"P"], function(m) rgamma(1000, shape=m/pars1["Pobs"], scale=pars1["Pobs"]))

par(mfrow=c(1,2), oma=rep(0.5,4), mar=c(5,5,0,0))
plot(out$time, out$Lobs, type='l', lwd=3, col="darkgrey", ylim=c(0,2), xlab="Age", ylab="Length (mm)", cex.lab=1.5)
for (i in 1:ncol(bootL)) lines(rep(unique(data$age)[i],2), sort(bootL[,i])[c(25,975)], lwd=3)
points(data$age, data$length, col=1, cex=1.5)

plot(out$time, out$P, type='l', lwd=3, col="darkgrey", ylim=c(0, 700000), xlab="Age", ylab="Parasite burden", cex.lab=1.5)
for (i in 1:ncol(bootP)) lines(rep(unique(data$age)[i],2), sort(bootP[,i])[c(25,975)], lwd=3)
points(data$age, data$spores, col=1, cex=1.5)


par(mfrow=c(1,2), oma=rep(0.5,4), mar=c(5,5,0,0))
plot(out1$time, out1$Lobs, type='l', lwd=3, col="darkgrey", ylim=c(0,2), xlab="Age", ylab="Length (mm)", cex.lab=1.5)
for (i in 1:ncol(bootL1)) lines(rep(unique(data$age)[i],2), sort(bootL1[,i])[c(25,975)], lwd=3)
points(data$age, data$length, col=1, cex=1.5)

plot(out1$time, out1$P, type='l', lwd=3, col="darkgrey", ylim=c(0, 700000), xlab="Age", ylab="Parasite burden", cex.lab=1.5)
for (i in 1:ncol(bootP1)) lines(rep(unique(data$age)[i],2), sort(bootP1[,i])[c(25,975)], lwd=3)
points(data$age, data$spores, col=1, cex=1.5)

