## Testing

source("Growth_reproduction_trajectory_matching_infecteds_model_2.R")

## Load the observed data
x <- read.csv("Cat_data/all_data.csv") %>% as.data.frame
x2 <- subset(x, infected==1)
data <- x2[,c("day.of.assay","length..mm.","Spore.per.animal")]
colnames(data) <- c("age","length","spores")
## how to deal with 0s?
data$spores[is.na(data$spores)] <- 0

## start w/a simple case: only estimate parasite parameters, assume all others equal to parasite free values
fixpars <- c(rho=0.152, K=0.68, km=0.073, v=10, F0=1e6/30, Lobs=0.0102)
estpars <- c(phimax=0.1, hP=1e3, eP=1e8, m=0.1, tau=1, Pobs=1000)

## unconstrained scale transformation
transform <- c('logit',rep('log',5))
parorder <- c("rho","K","km","v","F0","Lobs","phimax","hP","eP","m","tau","Pobs")

## combine into a single vector and sort
pars <- c(estpars, fixpars)
pars[match(parorder, names(pars))] -> pars
if(any(names(pars)!=parorder)) stop("parameter are in the wrong order")

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
        P=100,
        Z=0)
y0["E"] <- unname(y0["W"]*pars["rho"]/pars["v"])

## Simulate the system
try(ode(y0,
        times=seq(0,35,0.1),
        func="derivs",
        parms=pars,
        dllname="tm_deb_parasite_2",
        initfunc="initmod",
        events=list(data=eventdat))) -> out
## extract only the data points that can be compared against the true data
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
             dnorm(x=data$spores[data$age==d],
                   mean=pred$Z[pred$time==d],
                   sd=pars["Pobs"],
                   log=TRUE) %>% sum
             ) %>% sum
       ) %>% sum -> lik
print(lik)

## compare against trajectory matching code
tm_obj(estpars=par_transform(estpars, transform),
       data=data,
       fixpars=fixpars,
       parorder=parorder,
       transform=transform)

optimizer(estpars=estpars, fixpars=fixpars, parorder=parorder, transform=transform, obsdata=data, eval.only=TRUE, method="optim")

####################################################################
####################################################################
########    TRAJECTORY MATCHING OF THE REAL DATA!!!!    ############
####################################################################
####################################################################

source("Growth_reproduction_trajectory_matching_infecteds_model_2.R")

## parameters fixed at the uninfected data values
fixpars <- c(rho=0.152, K=0.68, km=0.073, v=10, F0=1e6/30, Lobs=0.0102)

## unconstrained scale transformation
transform <- c('logit',rep('log',5))
parorder <- c("rho","K","km","v","F0","Lobs","phimax","hP","eP","m","tau","Pobs")

estpars <- c(phimax=0.1, hP=1e3, eP=1e8, m=0.1, tau=1, Pobs=1000)

## Estimate the parameters
box <- cbind(lower=c(aP=0.01, hP=10, eP=1e5, m=1e-1, tau=0.1, Pobs=1e4),
             upper=c(aP=0.5, hP=1000, eP=1e10, m=10, tau=20, Pobs=1e6))
sobolDesign(lower=box[,'lower'],
            upper=box[,'upper'],
            nseq=300000) %>%
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
guesses[order(guess_lik)[1:3000]] -> refine
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
