## load all of the Daphnia data for genotype 3D at the highest food level
size.3d.aa = read.csv('~/Dropbox/Cressler/DEB_fitting/Adriana_Daphnia_data/3D-AA-GROWTH.csv', header=F, na.strings='na', colClasses='character')
birth.3d.aa = read.csv('~/Dropbox/Cressler/DEB_fitting/Adriana_Daphnia_data/3D-AA-ASEXUAL.csv', header=F, na.strings='na', colClasses='character')
## assign each individual to a position in a list
age <- as.numeric(size.3d.aa[6:nrow(size.3d.aa),1])
data <- vector(mode='list', length=(ncol(size.3d.aa)-1))
for (i in 1:length(data)) {
    print(i)
    Lobs <- as.numeric(size.3d.aa[6:nrow(size.3d.aa),i+1])
    end <- min(which(Lobs==-1))-1
    Lobs <- Lobs[1:end]
    Robs <- as.numeric(birth.3d.aa[6:nrow(birth.3d.aa),i+1])
    Robs <- Robs[1:end]
    data[[i]] = as.data.frame(cbind(age=age[1:end],cbind(Lobs,Robs)))
}

## use trajectory matching to find the best-fit parameter set for each
## individual, beginning from a global search
require(pomp)
require(plyr)

## create a single set of parameter guesses to try on each individual
## to begin the global search. Bill sent the parameter guesses he used
## to initialize his searches - I will use these as guides.
box <- cbind(
             lower=c(alpha=0.5,kappa=0.1,beta=0.01,
             W.scalar=5e-5,L.scalar=0.1,Rm.scalar=0.1,time.scalar=0.01,
             f=0.05,T_M=2,
             W.0=1,Lmeas.0=0.01,Rm.0=1e-9,L.0=0.01,
             L.sd=0.0001,R.sd=0.0001),
             upper=c(alpha=50,kappa=0.9,beta=5,
             W.scalar=5e-3,L.scalar=20,Rm.scalar=20,time.scalar=3,
             f=0.05,T_M=2,
             W.0=100,Lmeas.0=1,Rm.0=1e-9,L.0=1,
             L.sd=0.1,R.sd=0.1)
             )
set.seed(12340987)
guesses <- sobolDesign(lower=box[,'lower'],upper=box[,'upper'],nseq=10000)
guesses$id <- seq_len(nrow(guesses))
joblist <- dlply(guesses, ~id, unlist)
## strip out id
for (i in 1:length(joblist))
    joblist[[i]] <- joblist[[i]][-which(names(joblist[[i]])=='id')]

fix.pars <- c('f','T_M','Rm.0')
est.pars <- c('alpha','kappa','beta','W.scalar','L.scalar','Rm.scalar','time.scalar','W.0','Lmeas.0','L.0','L.sd','R.sd')
transform <- c('log','logit',rep('log',13)) ## rules for transforming to the estimation scale


source('trajectory.matching.R')
## Check to make sure everything appears to be working
this.data <- data[[1]]
this.parset <- joblist[[1]]
trajectory.matching(data=this.data, pars=this.parset, transform=transform,
                    est=est.pars, method='subplex',
                    odefunc="derivs", initfunc="initmod", dllname="nondim_deb",
                    eval.only=TRUE)

results = vector(mode='list', length=length(data))
for (i in 1:length(data)) {
    print(i)
    this.data <- data[[i]]
    if (nrow(this.data) > 7) {
        lik <- lapply(joblist, function(x) trajectory.matching(data=this.data,
                                                               pars=x,
                                                               transform=transform,
                                                               est=est.pars,
                                                               method='subplex',
                                                               odefunc="derivs",
                                                               initfunc="initmod",
                                                               dllname="nondim_deb",
                                                               eval.only=TRUE))
        this.results <- ldply(lik, function(x) c(x$par, loglik=x$lik))
        this.results <- this.results[order(this.results$loglik),]
        results[[i]] <- this.results
    }
}
save(results, file='~/Dropbox/Cressler/DEB_fitting/results-3D-AA-stage-one.rda')

## Go through a pull out only the top 25 parameter sets for each individual
good.params <- c()
for (i in 1:length(results)) {
    if (is.data.frame(results[[i]]))
        good.params <- rbind(good.params, results[[i]][1:25,])
}
## Look for duplicates by looking at .id (which is the number of the
## parameter set)
num.duplicates <- length(unique(good.params$.id[duplicated(good.params$.id)]))
## 218 parameter sets show up more than once

## Going back to the 10,000 parameter sets, count how many times each
## parameter set shows up in one of the top 25 fits.
n <- sapply(1:10000, function(x) sum(good.params$.id==x))
most.common.n <- rev(order(n))[1:10] ## all show up as one of the best fitting 25 parameter sets for >=19 individuals (of 58 individuals with enough data)
best.params <- c()
for (i in most.common.n)
    best.params <- rbind(best.params, joblist[[i]])

## Take the top 100 parameter sets for each individual and do a full
## trajectory.matching run, using subplex with only 500 steps initially.
results2 <- vector(mode='list', length=length(results))
for (i in 1:length(results)) {
    if (is.data.frame(results[[i]])) {
        this.data <- data[[i]]
        vals <- array(NA, dim=c(100, 18))
        for (j in 1:100) {
            print(j)
            this.parset <- results[[i]][j,2:16] ## extract the parameter values
            start.loglik <- results[[i]][j,17]
            x <- trajectory.matching(data=this.data,
                                     pars=this.parset,
                                     transform=transform,
                                     est=est.pars,
                                     method='subplex',
                                     odefunc="derivs",
                                     initfunc="initmod",
                                     dllname="nondim_deb",
                                     eval.only=FALSE,
                                     maxit=500)
            finish.loglik <- x$lik
            loglik.diff <- start.loglik-finish.loglik
            vals[j,] <- c(x$par, loglik=x$lik, conv=x$conv, dloglik=loglik.diff)
        }
        results2[[i]] <- vals
        save(results, file='~/Dropbox/Cressler/DEB_fitting/results-3D-AA-stage-two.rda')
    }
}


