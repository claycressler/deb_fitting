require(pomp)
require(plyr)

## Code for trajectory matching using the batch transfer ODE model
source('trajectory.matching.R')

## Generate a dataset with measurement error
true.pars <- c(kappa=0.6, km=0.33, eG=0.0017, eR=0.00868, nu=18.1, Rmbar=0.0189, pam=0.005187, Fh=3.09e-5, eA=0.7, vol=20, E.0=0.005187/18.1*0.85^3, L.0=0.85, Re.0=0, R.0=0, F.0=0.05, L.sd=0.02)
order <- c('kappa','km','eG','eR','nu','Rmbar','pam','Fh','eA','vol','E.0','L.0','Re.0','R.0','F.0','L.sd')
y0 <- c(y1=0.005187/18.1*0.85^3, y2=0.85, y3=0, y4=0, y5=0.05)
out <- ode(y0, times=seq(0,75), func="derivs", parms=true.pars, dllname="my_deb_batch_constantK", initfunc="initmod", events=list(func="event",time=seq(1,75,3)))
Lobs <- rnorm(length(out[(seq(1,75,3)+1),3]),out[(seq(1,75,3)+1),3],0.02)
Robs <- rpois(length(out[(seq(1,75,3)+1),5]),out[(seq(1,75,3)+1),5])
data <- as.data.frame(cbind(time=seq(1,75,3),cbind(Lobs,Robs)))

## Estimate the parameters of the fitted model
#x2 <- trajectory.matching(data=data,pars=allpars,est=c('kappa','km','nu'))

## An algorithm for performing a global search of parameter space.
## What parameters are we going to estimate?
estimated.pars <- c('kappa','km','eG','eR','nu','Rmbar','pam','Fh','eA','E.0','L.0')
## Encompass our initial ignorance of the estimated parameter values by a box
box <- cbind(
             lower=c(kappa=0.1, km=0.01, eG=0.0001, eR=0.0001, nu=1, Rmbar=0.001, pam=0.0005, Fh=8e-5, eA=0.5, E.0=0.00001, L.0=0.5),
             upper=c(kappa=0.9, km=1, eG=0.01, eR=0.01, nu=25, Rmbar=0.1, pam=0.05, Fh=8e-4, eA=1, E.0=0.0001, L.0=1.2)
             )
## Create a Sobol low discrepancy sequence of length 10000
guesses <- sobolDesign(lower=box[,'lower'],upper=box[,'upper'],nseq=10000)
guesses$id <- seq_len(nrow(guesses))
joblist <- dlply(guesses, ~id, unlist)
## fixed parameters will be constant
fixed.pars <- c(vol=20, Re.0=0, R.0=0, F.0=0.05, L.sd=0.02)
## calculate the logliklihood of each parameter combination
lik <- lapply(joblist, function(x) {
    allpars <- sort.params(x,fixed.pars,order);
    estpars <- names(x);
    estpars <- estpars[-which(estpars=='id')];
    trajectory.matching(data=data, pars=allpars, est=estpars, eval.only=TRUE)
})
## convert results to a dataframe
results = ldply(lik, function(x) c(x$par,loglik=x$lik))

save(data, file='Sample_batch_transfer_data.rda')
save(results, file='Sample_batch_data_stage_one.rda')

## sort by log-likelihood
nresults = results[order(results$loglik),]

## for each of the first 100 rows in nresults, run trajectory.matching
tm.results <- vector(mode='list', length=100)
for (i in 1:100) {
    y <- trajectory.matching(data=data,
                             pars=nresults[i,2:17],
                             est=estimated.pars,
                             method='Nelder-Mead',
                             maxit=1000)
    tm.results[[i]] <- y
}
tm.results <- ldply(tm.results, function(x) c(x$par, loglik=x$lik))
tm.results <- tm.results[order(tm.results$loglik),]
save(tm.results, file='Sample_batch_data_stage_two.rda')

## Looking at these results, it is interesting to note that the ones
## with the lowest negative log-likelihood are not necessarily those
## that were the best log-likelihoods previously.

## I am going to take the top 20, create novel variation around each
## putative optimal parameter set, and then run again, running
## Nelder-Mead for 2000 timesteps
tm.vary <- array(NA, dim=c(200,16))
for (q in 1:20) {
    tm.vary[((q-1)*10+1):(q*10),] <-
        t(sapply(1:10, function(x) {
            n <- names(tm.results[q,1:16])
            i <- which(n%in%estimated.pars)
            p <- as.numeric(tm.results[q,1:16])
            p[i] <- rnorm(length(p[i]),mean=p[i],sd=p[i]/10)
            p
        }))
}
colnames(tm.vary) <- colnames(tm.results)[1:16]
## Make sure all the parameters are positive
any(!(tm.vary >= 0))
## Double check to make sure that neither eA nor kappa is greater than
tm.vary[which(tm.vary[,'kappa'] > 1),'kappa'] <- 0.99
tm.vary[which(tm.vary[,'eA']>1),'eA'] <- 0.99

##
tm.vary.results <- vector(mode='list', length=200)
for (i in 1:200) {
    y <- try(trajectory.matching(data=data,
                             pars=tm.vary[i,],
                             est=estimated.pars,
                             method='Nelder-Mead',
                                 maxit=2000))
    if (inherits(y, 'try-error'))
        tm.vary.results[[i]] <- NA
    else
        tm.vary.results[[i]] <- y
    save(tm.vary.results, file='Sample_batch_data_stage_three_optim.rda')
}

tm.vary.results.2 <- vector(mode='list', length=200)
for (i in 1:200) {
    y <- try(trajectory.matching(data=data,
                                 pars=tm.vary[i,],
                                 est=estimated.pars,
                                 method='subplex'))
    if (inherits(y, 'try-error'))
        tm.vary.results.2[[i]] <- NA
    else
        tm.vary.results.2[[i]] <- y
    save(tm.vary.results.2, file='Sample_batch_data_stage_three_subplex.rda')
}

## compare the results from the two different fitting methods

## 26 of the parameter sets couldn't be fitted using
## trajectory.matching with Nelder-Mead, returning errors. Remove
## these.
err.ind <- which(sapply(tm.vary.results, function(x) is.na(unname(x[1]))))
good.ind <- seq(1,200)[-err.ind]
tm.vary.results.new <- vector(mode='list', length=200-length(err.ind))
for (i in 1:length(tm.vary.results.new))
    tm.vary.results.new[[i]] <- tm.vary.results[[good.ind[i]]]
tm.vary.results <- ldply(tm.vary.results.new, function(x) c(x$par, loglik=x$lik, conv=x$conv))
tm.vary.results.opt <- tm.vary.results[order(tm.vary.results$loglik),]
save(tm.vary.results.opt, file='Sample_batch_data_stage_three_optim.rda')

## 40 returned NA loglikelihoods - remove these
tm.vary.results.sub <- ldply(tm.vary.results.2, function(x) c(x$par, loglik=x$lik, conv=x$conv))
tm.vary.results.sub <- tm.vary.results.sub[order(tm.vary.results.sub$loglik),]
tm.vary.results.sub <- tm.vary.results.sub[1:159,]
save(tm.vary.results.sub, file='Sample_batch_data_stage_three_subplex.rda')


############################################################
## First off, compare the two best-fitting parameter sets
tm.vary.results.opt[1,]
tm.vary.results.sub[1,]
## These are totally, totally different. One has K=0.98, but is able
## to reproduce the data by taking eR to incredibly low values (so
## reproduction is easy) and increasing somatic maintenance and the
## cost of growth. Too many degrees of freedom gives the algorithm too
## much room to play.

## Hard to visualize unless I hack off some of the crappy fits.

## Focus on a few of the parameters: kappa, km, eG, eR, Rmbar, pam,
## eA, loglik. nu and Fh have lots of variation, making them hard to
## visualize, and I will leave off the two initial conditions

## Focus just on the subplex results for the moment
i <- which(names(tm.vary.results.sub)%in%estimated.pars)
pairs(tm.vary.results.sub[tm.vary.results.sub$loglik<45,i])

## Lots of covariation here, which suggests that the model could
## likely be reparamterized

## Which parameters are well estimated?
par(mfrow=c(3,4))
for (j in i) {
    hist(tm.vary.results.sub[tm.vary.results.sub$loglik<45,j],
         main=paste(names(tm.vary.results.sub)[j]),
         xlab='')
    abline(v=true.pars[names(tm.vary.results.sub)[j]], col=2)
}

quartz()
par(mfrow=c(3,4))
for (j in i) {
    hist(tm.vary.results.opt[tm.vary.results.opt$loglik<45,j],
         main=paste(names(tm.vary.results.opt)[j]),
         xlab='')
    abline(v=true.pars[names(tm.vary.results.opt)[j]], col=2)
}



## Are combinations of parameters well-estimated?
## Try Rmbar/eR as a combination:
hist(tm.vary.results.sub[tm.vary.results.sub$loglik<45,'Rmbar']/tm.vary.results.sub[tm.vary.results.sub$loglik<45,'eR'], main=expression('subplex Estimates'~~'of'~~bar(R[m])/epsilon[R]),xlab='')
abline(v=true.pars['Rmbar']/true.pars['eR'],col=2)
text(x=(true.pars['Rmbar']/true.pars['eR']), y=0.5, pos=1, labels='Truth', col=2)

## Try kappa/eG
hist(tm.vary.results.sub[tm.vary.results.sub$loglik<45,'kappa']/tm.vary.results.sub[tm.vary.results.sub$loglik<45,'eG'], main=expression('subplex Estimates'~~'of'~~kappa/epsilon[G]),xlab='')
abline(v=true.pars['kappa']/true.pars['eG'],col=2)
text(x=(true.pars['kappa']/true.pars['eG']), y=0.5, pos=1, labels='Truth', col=2)

## Try kappa/eR
hist(tm.vary.results.sub[tm.vary.results.sub$loglik<45,'kappa']/tm.vary.results.sub[tm.vary.results.sub$loglik<45,'eR'], main=expression('subplex Estimates'~~'of'~~kappa/epsilon[R]),xlab='',breaks=100000,xlim=c(0,100))
abline(v=true.pars['kappa']/true.pars['eR'],col=2)
text(x=(true.pars['kappa']/true.pars['eR']), y=0.5, pos=1, labels='Truth', col=2)

attach(tm.vary.results.sub)
estnames <- names(tm.vary.results.sub[,i])
iseq <- seq(1,length(estnames))
par(mfrow=c(length(i),length(i)), mar=rep(1,4), oma=rep(0,4))
for (j in iseq) {
    for (k in iseq) {
        if (j==k) {
            plot(1,type='n',axes=F,xlab='',ylab='')
            text(x=0.5,y=0.5,labels=estnames[j])
        }
        else {
            hist(tm.vary.results.sub[loglik<45,i[j]]/tm.vary.results.sub[loglik<45,i[k]], xlab='', main=paste0(estnames[j],'/',estnames[k]))
            abline(v=true.pars[estnames[j]]/true.pars[estnames[k]], col=2)
            #text(x=true.pars[estnames[j]]/true.pars[estnames[k]],y=0.5,labels='Truth')
        }
    }
}

