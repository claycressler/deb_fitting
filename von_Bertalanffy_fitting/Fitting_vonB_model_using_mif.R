## The likelihood of observing length Vobs and births Robs given that
## L and R were the actual length and births
" lik = dnorm(Lobs,L,L_sd,give_log); " -> dmeas

## Simulate a measurement
" Lobs = rnorm(L,L_sd); " -> rmeas

## Encode the model process simulator
## G_sd is stochastic variation in growth due to, e.g., stochastic
## variation in food supply
' double rate; //transition rates
  double trans; // transition numbers
  rate = r*(Linf-L);
  trans = rnorm(rate*dt, G_sd*sqrt(dt));
  L += trans;
' -> stepfn

' double rate; //transition rates
  rate = r*(Linf-L);
  DL += rate;
' -> skel

"
  Tr = exp(r);
  TLinf = exp(Linf);
  TL_sd = exp(L_sd);
  TL_0 = exp(L_0);
  TG_sd = exp(G_sd);
" -> trans

"
  Tr = log(r);
  TLinf = log(Linf);
  TL_sd = log(L_sd);
  TL_0 = log(L_0);
  TG_sd = log(G_sd);
" -> untrans

## Load the data for Glen genotype
load('~/data/Research/DEB_fitting/Adriana_Daphnia_data/size.glen.rda')
size.glen = size.glen[,1:(ncol(size.glen)-1)]

####################### EXPLORATION ##############################
## For super-simplicity, just look at the very first set of data
trial.glen = size.glen[,1:2]
colnames(trial.glen) = c('Age','Lobs') ## Lobs is name in DEB_growth.R

pompBuilder(name='VB',
            data=as.data.frame(trial.glen),
            times='Age',
            t0=0,
            dmeasure=dmeas,
            rmeasure=rmeas,
            step.fn=stepfn,
            step.fn.delta.t=1/24,
            skeleton.type='vectorfield',
            skeleton=skel,
            parameter.transform=trans,
            parameter.inv.transform=untrans,
            statenames=c('L'),
            paramnames=c('r','Linf','L.0','L.sd','G.sd')
            ) -> VB

coef(VB) <- c(r=0.05, Linf=2.75, L.0=0.75, L.sd=0.04, G.sd=0.04)
x <- simulate(VB, states=T, as.data.frame=T)
plot(Lobs~Age, data=trial.glen, type='l')
lines(L~time, data=x, type='l', col=2)

pf1 <- pfilter(VB,Np=1000)

mf1 <- mif(VB,
           Nmif=200,
           transform=T,
           var.factor=2,
           cooling.factor=0.99,
           ic.lag=25,
           Np=1000,
           ivps=c('L.0'),
           rw.sd=c(L.0=0.05,r=0.02,Linf=0.02,L.sd=0.02,G.sd=0.02)
           )

pf2 <- pfilter(mf1, Np=1000)


#########################################################################

## Develop a strategy for global searches, initializing at many random
## starts, refining, etc. that I feel confident gives me a global
## maximum.
require(plyr)
## Encompass my ignorance of parameter values in a box
box <- cbind(
             lower=c(r=0.001,Linf=2,L.0=0.25,L.sd=0.001,G.sd=0.001),
             upper=c(r=1,Linf=3,L.0=1,L.sd=1,G.sd=1)
             )
## Generate many points well-distributed through the box using Sobol
## low-discrepacy sequences
guesses <- sobolDesign(lower=box[,'lower'],upper=box[,'upper'],nseq=1000)
guesses$id <- seq_len(nrow(guesses))

## Turn the guesses into a list that can be passed to lapply
joblist <- dlply(guesses,~id,unlist)

workfn <- function (start, po, Np, ...) {
    ## particle filter to calculate log-likelihood for each of the
    ## parameter sets in guesses
    pf <- try(pfilter(po,params=start,Np=Np,...))
    if (inherits(pf,"try-error")) {
        ll <- NA
        params <- start
    } else {
        ll <- logLik(pf)
        params <- coef(pf)
    }
    list(
         start=start,
         params=params,
         loglik=ll
         )
}

## run the particle filter on the pomp objects VB1-VB9
for (i in 1:9) {
    assign(paste('results',i,sep=''),lapply(joblist,workfn,po=get(paste('VB',i,sep='')),Np=1000))
}

results <- lapply(joblist,workfn,
## So create a bunch of mif objects fitting each of the datasets
## independently for a bit, then start threading between them using
## continue and substituting the shared parameter values among the mif
## objects.

for (i in 2:ncol(size.glen)) {
    dat <- as.data.frame(size.glen[,c(1,i)])
    ## Remove rows with missing datapoints (NAs), as mif will produce an error.
    if (any(is.na(dat)))
        dat <- dat[-unname(which(is.na(dat),arr.ind=T)[1,1]),]
    colnames(dat) <- c('Age','Lobs')
    assign(paste('VB',i-1,sep=''),
           pompBuilder(name=paste('VB',i-1,sep=''),
                       data=dat,
                       times='Age',
                       t0=0,
                       dmeasure=dmeas,
                       rmeasure=rmeas,
                       step.fn=stepfn,
                       step.fn.delta.t=1/24,
                       skeleton.type='vectorfield',
                       skeleton=skel,
                       parameter.transform=trans,
                       parameter.inv.transform=untrans,
                       statenames=c('L'),
                       paramnames=c('r','Linf','L.0','L.sd','G.sd')
                       )
           )
    assign(paste('mif',i-1,sep=''),
           mif(get(paste('VB',i-1,sep='')),
               Nmif=1,
               start=c(r=0.05,Linf=2.75,L.0=0.75,L.sd=0.04,G.sd=0.04),
               transform=T,
               var.factor=2,
               cooling.factor=0.99,
               ic.lag=nrow(dat),
               Np=1000,
               ivps=c('L.0'),
               rw.sd=c(L.0=0.05,r=0.02,Linf=0.02,L.sd=0.02,G.sd=0.02)
               )
           )
}

## Now, continue for Nmif times, threading between the different datasets - all possible combinations
Nmif <- 5
for (m in 1:Nmif) {
    print(m)
    shared.params <- c('r','Linf','L.0','L.sd','G.sd')
    mif.objects <- paste('mif',seq(1,ncol(size.glen)-1),sep='')
    ## trade coefficients sets between all possible combinations
    trade.order <- expand.grid(mif.objects,mif.objects,stringsAsFactors=F)
    ## remove rows where trade partners are identical
    trade.order <- trade.order[-which(apply(trade.order,1,function(x)x[1]==x[2])),]
    ## randomize the order of substitutions each mif iteration
    trade.order <- trade.order[sample(1:nrow(trade.order),nrow(trade.order),replace=F),]

    ## continue() each mif, trading parameter values among the different datasets
    for (row in 1:nrow(trade.order)) {
        mif.obj1 <- get(unname(trade.order[row,1]))
        mif.obj2 <- get(unname(trade.order[row,2]))
        ## trade coefficients
        coef(mif.obj1,shared.params) <- coef(mif.obj2,shared.params)
        ## continue mif.obj1
        mif.obj1 <- continue(mif.obj1, Nmif=1)
        ## trade coefficients back
        coef(mif.obj2,shared.params) <- coef(mif.obj1,shared.params)
        ## continue mif.obj2
        assign(unname(trade.order[row,2]), continue(mif.obj2, Nmif=1))
        ## reassign
        assign(unname(trade.order[row,1]),mif.obj1)
        assign(unname(trade.order[row,2]),mif.obj2)
    }
}



## How to deal with multiple datasets.
## imagine having two pomp objects, po1 and po2, that have shared parameter values between them, but have independent dynamics. (If they have shared dynamics, you have to vectorize the state space instead.)
for (m in 1:Nmif) {
    shared.params <- c('shared parameters between po1, po2')
    coef(po1, shared.params) <- coef(po2, shared.params) ## Trade the shared parameter values between po1 and po2
    po1 <- continue(po1, Nmif=1)
    coef(po2, shared.params) <- coef(po1, shared.params) ## Trade the shared parameter values between po1 and po2
    po2 <- continue(po2, Nmif=1)
}

