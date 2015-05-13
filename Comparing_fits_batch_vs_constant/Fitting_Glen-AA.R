## load all of the Daphnia data for genotype GLEN at the highest food level
size.GLEN.aa = read.csv('~/Dropbox/Cressler/DEB_fitting/Adriana_Daphnia_data/GLEN-AA-GROWTH.csv', header=F, na.strings='na', colClasses='character')
birth.GLEN.aa = read.csv('~/Dropbox/Cressler/DEB_fitting/Adriana_Daphnia_data/GLEN-AA-ASEXUAL.csv', header=F, na.strings='na', colClasses='character')
## assign each individual to a position in a list
age <- as.numeric(size.GLEN.aa[6:nrow(size.GLEN.aa),1])
data <- vector(mode='list', length=(ncol(size.GLEN.aa)-1))
for (i in 1:length(data)) {
    Lobs <- as.numeric(size.GLEN.aa[6:nrow(size.GLEN.aa),i+1])
    end <- min(which(Lobs==-1))-1
    Lobs <- Lobs[1:end]
    Robs <- as.numeric(birth.GLEN.aa[6:nrow(birth.GLEN.aa),i+1])
    Robs <- Robs[1:end]
    data[[i]] = as.data.frame(cbind(age=age[1:end],cbind(Lobs,Robs)))
}
save(data, file='GLEN-AA_all_growth_egg_data.rda')

## use trajectory matching to find the best-fit parameter set for each
## individual, beginning from a global search
require(pomp)
require(plyr)
source('trajectory.matching.R')
name.order <- c('K','km','eG','eR','v','Rmbar','Imax','Fh','eA','E.0','L.0','R.0','F.0','L.sd')

#for (mod in c("batch", "constant")) {
mod <- 'constant'
    ## Step 1: Evaluate many randomly chosen parameter sets
    ## Parameters to estimate:
    estimated.pars <- c('K','km','eG','eR','v','Rmbar','Imax','Fh','eA','E.0','L.0')
    ## Bound our ignorance of these parameters with a box
    box <- cbind(
                 lower=c(K=0.1, km=0.01, eG=0.0001, eR=0.0001, v=1, Rmbar=0.001, Imax=0.0005, Fh=0.01, eA=0.5, E.0=0.00001, L.0=0.5),
                 upper=c(K=0.9, km=1, eG=0.01, eR=0.01, v=25, Rmbar=0.1, Imax=0.05, Fh=0.5, eA=1, E.0=0.0001, L.0=1.2)
                 )
    set.seed(12340987)
    ## Create 1000 random guesses at the parameters
    guesses <- sobolDesign(lower=box[,'lower'],upper=box[,'upper'],nseq=1000)
    guesses$id <- seq_len(nrow(guesses))
    joblist <- dlply(guesses, ~id, unlist)
    ## strip out id
    for (i in 1:length(joblist))
        joblist[[i]] <- joblist[[i]][-which(names(joblist[[i]])=='id')]

    ## Parameters whose value is fixed (note: F is measured in mgC/L):
    fixed.pars <- c(R.0=0,F.0=0.1/0.02,L.sd=0.02)

    ## Calculate the likelihood of each parameter set
    results = vector(mode='list', length=length(data))
    for (i in 1:length(data)) {
        this.data <- data[[i]]
        ## only calculate likelihoods for animals that lived at least 15 days
        if (nrow(this.data) > 7) {
            lik <- lapply(joblist, function(x) trajectory.matching(data=this.data,
                                                                   fixed.pars=fixed.pars,
                                                                   estd.pars=x,
                                                                   method='subplex',
                                                                   eval.only=TRUE,
                                                                   model=mod))
            this.results <- ldply(lik, function(x) c(x$par, loglik=x$lik))
            this.results <- this.results[order(this.results$loglik),]
            results[[i]] <- this.results
        }
    }
    save(results, file=paste0('~/Dropbox/Cressler/DEB_fitting/Comparing_fits_batch_vs_constant/',mod,'_results-GLEN-AA-stage-one.rda'))

    ## Step 2:
    ## Take the top 100 parameter sets for each individual and do a full
    ## trajectory.matching run, using subplex with only 500 steps initially.
    results2 <- vector(mode='list', length=length(results))
    for (i in 31:length(results)) {
        if (is.data.frame(results[[i]])) {
            this.data <- data[[i]]
            vals <- array(NA, dim=c(100, 17))
            for (j in 1:100) {
                print(paste('dataset',i,'parameter set',j))
                this.parset <- results[[i]][j,2:15] ## extract the parameter values
                estd.pars <- this.parset[-which(names(this.parset)%in%names(fixed.pars))]
                start.loglik <- results[[i]][j,16]
                x <- trajectory.matching(data=this.data,
                                         fixed.pars=fixed.pars,
                                         estd.pars=estd.pars,
                                         method='subplex',
                                         eval.only=FALSE,
                                         model=mod,
                                         maxit=500)
                finish.loglik <- x$lik
                loglik.diff <- start.loglik-finish.loglik
                vals[j,] <- c(x$par, loglik=x$lik, conv=x$conv, dloglik=loglik.diff)
            }
            vals <- as.data.frame(vals)
            colnames(vals) <- c(name.order,'loglik','conv','dloglik')
            results2[[i]] <- vals
            save(results2, file=paste0('~/Dropbox/Cressler/DEB_fitting/Comparing_fits_batch_vs_constant/',mod,'results-GLEN-AA-stage-two.rda'))
        }
    }



## Step 3:
## Further refine
#    results3 <- vector(mode='list', length=length(results2))
#    for (i in 1:length(results2)) {
#        if (is.data.frame(results2[[i]])) {
#            this.data <- data[[i]]
#            vals <- array(NA, dim=c(100, 18))
#            for (j in 1:100) {
#                if (!is.na(results2[[i]]$loglik[j])) {
#                    print(paste('dataset',i,'parset',j))
#                    this.parset <- results2[[i]][j,1:15] ## extract the parameter values
#                    start.loglik <- results2[[i]]$loglik[j]
#                    x <- try(trajectory.matching(data=this.data,
#                                                 pars=this.parset,
#                                                 transform=transform,
#                                                 est=est.pars,
#                                                 method='subplex',
#                                                 odefunc="derivs",
#                                                 initfunc="initmod",
#                                                 dllname="nondim_deb",
#                                                 eval.only=FALSE,
#                                                 maxit=1500))
#                    if (inherits(x,'try-error'))
###                        vals[j,] <- c(as.numeric(this.parset), loglik=start.loglik, conv=NA, dloglik=0)
#                    else {
#                        finish.loglik <- x$lik
#                        loglik.diff <- start.loglik-finish.loglik
#                        vals[j,] <- c(x$par, loglik=x$lik, conv=x$conv, dloglik=loglik.diff)
#                    }
#                }
#                else vals[j,] <- as.numeric(results2[[i]][j,])
#            }
#            vals <- as.data.frame(vals)
#            colnames(vals) <- c('alpha','kappa','beta','W.scalar','L.scalar','Rm.scalar','time.scalar','f','T_M','W.0','Lmeas.0','Rm.0','L.0','L.sd','R.sd','loglik','conv','dloglik')
#            results3[[i]] <- vals
#            save(results3, file='~/Dropbox/Cressler/DEB_fitting/results-GLEN-AA-stage-three.rda')
#        }
#    }
#}

