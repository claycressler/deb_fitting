dyn.load('Bill_deb_redimensionalized.so')
require(deSolve)

## 2C-AA ##
pars.2c <- read.csv('2C-AA-p.firstrun-food_0.05.csv')
## the first row appears to be garbage, based on comparison with the
## real data
pars.2c <- pars.2c[-1,]
## load the real data
load('2C-AA_all_growth_egg_data.rda')
## add age at death to each parameter vector
pars.2c$aod <- sapply(data, function(x) x[nrow(x),'age'])
## double-check to make sure everything works out - all rows with NA
## parameter values should have age at death less than 16
cbind(pars.2c$alpha, pars.2c$aod)
## remove all NA rows
pars.2c <- pars.2c[-which(is.na(pars.2c[,1])),]
## transform the parameters to the natural scale
pars.2c[,'kappa'] <- exp(pars.2c[,'kappa'])/(1+exp(pars.2c[,'kappa']))
pars.2c[,-which(names(pars.2c)=='kappa')] <- exp(pars.2c[,-which(names(pars.2c)=='kappa')])
pars.2c[,'aod'] <- log(pars.2c[,'aod'])
## parameters needed by Bill_deb_redimensionalized
sim.pars.2c <- cbind(pars.2c[,c('alpha','kappa','beta')],
                     rep(0.1,nrow(pars.2c)),
                     rep(2,nrow(pars.2c)))
colnames(sim.pars.2c)[4:5] <- c('f','T_M')
x.2c.all <- vector(mode='list', length=nrow(pars.2c))
for (i in 1:nrow(pars.2c)) {
    ## length of time to run the simulation
    sim.times.2c <- seq(0,pars.2c[i,'aod'],0.1)*pars.2c[i,'d']
    ## initial conditions
    sim.y.2c <- c(pars.2c[i,'Wo.dim']*pars.2c[i,'d']/pars.2c[i,'a'],
                  pars.2c[i,'Lo.dim']*pars.2c[i,'d']/pars.2c[i,'b'],
                  0,
                  pars.2c[i,'Lo.dim']*pars.2c[i,'d']/pars.2c[i,'b'])
    names(sim.y.2c) <- c('W','Lmeas','R','L')
    ## simulate the system
    x.2c <- ode(y=sim.y.2c, times=sim.times.2c, func='derivs',
                parms=as.numeric(sim.pars.2c[i,]),
                dllname='Bill_deb_redimensionalized', initfunc='initmod')
    ## redimensionalize the output
    x.2c <- as.data.frame(x.2c)
    x.2c$time <- x.2c$time/pars.2c[i,'d']
    x.2c$Lmeas <- x.2c$Lmeas*pars.2c[i,'b']/pars.2c[i,'d']
    x.2c$R <- x.2c$R-pars.2c[i,'beta']
    x.2c$R[which(x.2c$R<=0)] <- 0
    x.2c$R <- x.2c$R*pars.2c[i,'c']/pars.2c[i,'d']
    x.2c.all[[i]] <- x.2c
}

## 4B-AA ##
pars.4b <- read.csv('4B-AA-p.firstrun-food_0.05.csv')
## the first row appears to be garbage, based on comparison with the
## real data
pars.4b <- pars.4b[-1,]
## load the real data
load('4B-AA_all_growth_egg_data.rda')
## add age at death to each parameter vector
pars.4b$aod <- sapply(data, function(x) x[nrow(x),'age'])
## double-check to make sure everything works out - all rows with NA
## parameter values should have age at death less than 16
cbind(pars.4b$alpha, pars.4b$aod)
## remove all NA rows
pars.4b <- pars.4b[-which(is.na(pars.4b[,1])),]
## transform the parameters to the natural scale
pars.4b[,'kappa'] <- exp(pars.4b[,'kappa'])/(1+exp(pars.4b[,'kappa']))
pars.4b[,-which(names(pars.4b)=='kappa')] <- exp(pars.4b[,-which(names(pars.4b)=='kappa')])
pars.4b[,'aod'] <- log(pars.4b[,'aod'])
## parameters needed by Bill_deb_redimensionalized
sim.pars.4b <- cbind(pars.4b[,c('alpha','kappa','beta')],
                     rep(0.1,nrow(pars.4b)),
                     rep(2,nrow(pars.4b)))
colnames(sim.pars.4b)[4:5] <- c('f','T_M')
x.4b.all <- vector(mode='list', length=nrow(pars.4b))
for (i in 1:nrow(pars.4b)) {
    ## length of time to run the simulation
    sim.times.4b <- seq(0,pars.4b[i,'aod'],0.1)*pars.4b[i,'d']
    ## initial conditions
    sim.y.4b <- c(pars.4b[i,'Wo.dim']*pars.4b[i,'d']/pars.4b[i,'a'],
                  pars.4b[i,'Lo.dim']*pars.4b[i,'d']/pars.4b[i,'b'],
                  0,
                  pars.4b[i,'Lo.dim']*pars.4b[i,'d']/pars.4b[i,'b'])
    names(sim.y.4b) <- c('W','Lmeas','R','L')
    ## simulate the system
    x.4b <- ode(y=sim.y.4b, times=sim.times.4b, func='derivs',
                parms=as.numeric(sim.pars.4b[i,]),
                dllname='Bill_deb_redimensionalized', initfunc='initmod')
    ## redimensionalize the output
    x.4b <- as.data.frame(x.4b)
    x.4b$time <- x.4b$time/pars.4b[i,'d']
    x.4b$Lmeas <- x.4b$Lmeas*pars.4b[i,'b']/pars.4b[i,'d']
    x.4b$R <- x.4b$R-pars.4b[i,'beta']
    x.4b$R[which(x.4b$R<=0)] <- 0
    x.4b$R <- x.4b$R*pars.4b[i,'c']/pars.4b[i,'d']
    x.4b.all[[i]] <- x.4b
}

## GLEN-AA ##
pars.GLEN <- read.csv('GLEN-AA-p.firstrun-food_0.05.csv')
## the first row appears to be garbage, based on comparison with the
## real data
pars.GLEN <- pars.GLEN[-1,]
## load the real data
load('GLEN-AA_all_growth_egg_data.rda')
## add age at death to each parameter vector
pars.GLEN$aod <- sapply(data, function(x) x[nrow(x),'age'])
## double-check to make sure everything works out - all rows with NA
## parameter values should have age at death less than 16
cbind(pars.GLEN$alpha, pars.GLEN$aod)
## remove all NA rows
pars.GLEN <- pars.GLEN[-which(is.na(pars.GLEN[,1])),]
## transform the parameters to the natural scale
pars.GLEN[,'kappa'] <- exp(pars.GLEN[,'kappa'])/(1+exp(pars.GLEN[,'kappa']))
pars.GLEN[,-which(names(pars.GLEN)=='kappa')] <- exp(pars.GLEN[,-which(names(pars.GLEN)=='kappa')])
pars.GLEN[,'aod'] <- log(pars.GLEN[,'aod'])
## parameters needed by Bill_deb_redimensionalized
sim.pars.GLEN <- cbind(pars.GLEN[,c('alpha','kappa','beta')],
                     rep(0.1,nrow(pars.GLEN)),
                     rep(2,nrow(pars.GLEN)))
colnames(sim.pars.GLEN)[4:5] <- c('f','T_M')
x.GLEN.all <- vector(mode='list', length=nrow(pars.GLEN))
for (i in 1:nrow(pars.GLEN)) {
    ## length of time to run the simulation
    sim.times.GLEN <- seq(0,pars.GLEN[i,'aod'],0.1)*pars.GLEN[i,'d']
    ## initial conditions
    sim.y.GLEN <- c(pars.GLEN[i,'Wo.dim']*pars.GLEN[i,'d']/pars.GLEN[i,'a'],
                  pars.GLEN[i,'Lo.dim']*pars.GLEN[i,'d']/pars.GLEN[i,'b'],
                  0,
                  pars.GLEN[i,'Lo.dim']*pars.GLEN[i,'d']/pars.GLEN[i,'b'])
    names(sim.y.GLEN) <- c('W','Lmeas','R','L')
    ## simulate the system
    x.GLEN <- ode(y=sim.y.GLEN, times=sim.times.GLEN, func='derivs',
                parms=as.numeric(sim.pars.GLEN[i,]),
                dllname='Bill_deb_redimensionalized', initfunc='initmod')
    ## redimensionalize the output
    x.GLEN <- as.data.frame(x.GLEN)
    x.GLEN$time <- x.GLEN$time/pars.GLEN[i,'d']
    x.GLEN$Lmeas <- x.GLEN$Lmeas*pars.GLEN[i,'b']/pars.GLEN[i,'d']
    x.GLEN$R <- x.GLEN$R-pars.GLEN[i,'beta']
    x.GLEN$R[which(x.GLEN$R<=0)] <- 0
    x.GLEN$R <- x.GLEN$R*pars.GLEN[i,'c']/pars.GLEN[i,'d']
    x.GLEN.all[[i]] <- x.GLEN
}

pars.RAMSEY <- read.csv('RAMSEY-AA-p.firstrun-food_0.05.csv')
## the first row appears to be garbage, based on comparison with the
## real data
pars.RAMSEY <- pars.RAMSEY[-1,]
## load the real data
load('RAMSEY-AA_all_growth_egg_data.rda')
## add age at death to each parameter vector
pars.RAMSEY$aod <- sapply(data, function(x) x[nrow(x),'age'])
## double-check to make sure everything works out - all rows with NA
## parameter values should have age at death less than 16
cbind(pars.RAMSEY$alpha, pars.RAMSEY$aod)
## remove all NA rows
pars.RAMSEY <- pars.RAMSEY[-which(is.na(pars.RAMSEY[,1])),]
## transform the parameters to the natural scale
pars.RAMSEY[,'kappa'] <- exp(pars.RAMSEY[,'kappa'])/(1+exp(pars.RAMSEY[,'kappa']))
pars.RAMSEY[,-which(names(pars.RAMSEY)=='kappa')] <- exp(pars.RAMSEY[,-which(names(pars.RAMSEY)=='kappa')])
pars.RAMSEY[,'aod'] <- log(pars.RAMSEY[,'aod'])
## parameters needed by Bill_deb_redimensionalized
sim.pars.RAMSEY <- cbind(pars.RAMSEY[,c('alpha','kappa','beta')],
                     rep(0.1,nrow(pars.RAMSEY)),
                     rep(2,nrow(pars.RAMSEY)))
colnames(sim.pars.RAMSEY)[4:5] <- c('f','T_M')
x.RAMSEY.all <- vector(mode='list', length=nrow(pars.RAMSEY))
for (i in 1:nrow(pars.RAMSEY)) {
    ## length of time to run the simulation
    sim.times.RAMSEY <- seq(0,pars.RAMSEY[i,'aod'],0.1)*pars.RAMSEY[i,'d']
    ## initial conditions
    sim.y.RAMSEY <- c(pars.RAMSEY[i,'Wo.dim']*pars.RAMSEY[i,'d']/pars.RAMSEY[i,'a'],
                  pars.RAMSEY[i,'Lo.dim']*pars.RAMSEY[i,'d']/pars.RAMSEY[i,'b'],
                  0,
                  pars.RAMSEY[i,'Lo.dim']*pars.RAMSEY[i,'d']/pars.RAMSEY[i,'b'])
    names(sim.y.RAMSEY) <- c('W','Lmeas','R','L')
    ## simulate the system
    x.RAMSEY <- ode(y=sim.y.RAMSEY, times=sim.times.RAMSEY, func='derivs',
                parms=as.numeric(sim.pars.RAMSEY[i,]),
                dllname='Bill_deb_redimensionalized', initfunc='initmod')
    ## redimensionalize the output
    x.RAMSEY <- as.data.frame(x.RAMSEY)
    x.RAMSEY$time <- x.RAMSEY$time/pars.RAMSEY[i,'d']
    x.RAMSEY$Lmeas <- x.RAMSEY$Lmeas*pars.RAMSEY[i,'b']/pars.RAMSEY[i,'d']
    x.RAMSEY$R <- x.RAMSEY$R-pars.RAMSEY[i,'beta']
    x.RAMSEY$R[which(x.RAMSEY$R<=0)] <- 0
    x.RAMSEY$R <- x.RAMSEY$R*pars.RAMSEY[i,'c']/pars.RAMSEY[i,'d']
    x.RAMSEY.all[[i]] <- x.RAMSEY
}

## Plots ##
par(mfrow=c(1,2), oma=rep(0.5,4), mar=c(5,5,1,1))
plot.new()
plot.window(xlim=c(0,90), ylim=c(0.5,3))
axis(1)
axis(2)
box('plot')
mtext(side=1, line=3, 'Age')
mtext(side=2, line=3, 'Length')
for (i in 1:length(x.2c.all))
    points(x.2c.all[[i]]$time, x.2c.all[[i]]$Lmeas, type='l', lwd=1.5, col='blue')
for (i in 1:length(x.4b.all))
    points(x.4b.all[[i]]$time, x.4b.all[[i]]$Lmeas, type='l', lwd=1.5, col='green')
for (i in 1:length(x.GLEN.all))
    points(x.GLEN.all[[i]]$time, x.GLEN.all[[i]]$Lmeas, type='l', lwd=1.5, col='gold')
for (i in 1:length(x.RAMSEY.all))
    points(x.RAMSEY.all[[i]]$time, x.RAMSEY.all[[i]]$Lmeas, type='l', lwd=1.5, col='red')

plot.new()
plot.window(xlim=c(0,90), ylim=c(0,60))
axis(1)
axis(2)
box('plot')
mtext(side=1, line=3, 'Age')
mtext(side=2, line=3, 'Cum. Eggs')
for (i in 1:length(x.2c.all))
    points(x.2c.all[[i]]$time, x.2c.all[[i]]$R, type='l', lwd=1.5, col='blue')
for (i in 1:length(x.4b.all))
    points(x.4b.all[[i]]$time, x.4b.all[[i]]$R, type='l', lwd=1.5, col='green')
for (i in 1:length(x.GLEN.all))
    points(x.GLEN.all[[i]]$time, x.GLEN.all[[i]]$R, type='l', lwd=1.5, col='gold')
for (i in 1:length(x.RAMSEY.all))
    points(x.RAMSEY.all[[i]]$time, x.RAMSEY.all[[i]]$R, type='l', lwd=1.5, col='red')

dev.copy2pdf(file='Bill_all_fit_trajectories.pdf')

## Calculate the mean trajectories for each genotype
## Determine the length of the longest-lived individual
nrows <- max(sapply(x.2c.all, function(x) nrow(x)))
ncols <- length(x.2c.all)
mean.2c.age <- seq(0,(nrows-1)*0.1,0.1)
all.2c.Lmeas <- array(NA, dim=c(nrows,ncols))
all.2c.R <- array(NA, dim=c(nrows,ncols))
for (i in 1:ncols) {
    this.Lmeas <- x.2c.all[[i]]$Lmeas
    all.2c.Lmeas[1:length(this.Lmeas),i] <- this.Lmeas
    this.R <- x.2c.all[[i]]$R
    all.2c.R[1:length(this.R),i] <- this.R
}
mean.2c.Lmeas <- apply(all.2c.Lmeas, 1, function(x) mean(x, na.rm=T))
mean.2c.R <- apply(all.2c.R, 1, function(x) mean(x, na.rm=T))

nrows <- max(sapply(x.4b.all, function(x) nrow(x)))
ncols <- length(x.4b.all)
mean.4b.age <- seq(0,(nrows-1)*0.1,0.1)
all.4b.Lmeas <- array(NA, dim=c(nrows,ncols))
all.4b.R <- array(NA, dim=c(nrows,ncols))
for (i in 1:ncols) {
    this.Lmeas <- x.4b.all[[i]]$Lmeas
    all.4b.Lmeas[1:length(this.Lmeas),i] <- this.Lmeas
    this.R <- x.4b.all[[i]]$R
    all.4b.R[1:length(this.R),i] <- this.R
}
mean.4b.Lmeas <- apply(all.4b.Lmeas, 1, function(x) mean(x, na.rm=T))
mean.4b.R <- apply(all.4b.R, 1, function(x) mean(x, na.rm=T))

nrows <- max(sapply(x.GLEN.all, function(x) nrow(x)))
ncols <- length(x.GLEN.all)
mean.GLEN.age <- seq(0,(nrows-1)*0.1,0.1)
all.GLEN.Lmeas <- array(NA, dim=c(nrows,ncols))
all.GLEN.R <- array(NA, dim=c(nrows,ncols))
for (i in 1:ncols) {
    this.Lmeas <- x.GLEN.all[[i]]$Lmeas
    all.GLEN.Lmeas[1:length(this.Lmeas),i] <- this.Lmeas
    this.R <- x.GLEN.all[[i]]$R
    all.GLEN.R[1:length(this.R),i] <- this.R
}
mean.GLEN.Lmeas <- apply(all.GLEN.Lmeas, 1, function(x) mean(x, na.rm=T))
mean.GLEN.R <- apply(all.GLEN.R, 1, function(x) mean(x, na.rm=T))

nrows <- max(sapply(x.RAMSEY.all, function(x) nrow(x)))
ncols <- length(x.RAMSEY.all)
mean.RAMSEY.age <- seq(0,(nrows-1)*0.1,0.1)
all.RAMSEY.Lmeas <- array(NA, dim=c(nrows,ncols))
all.RAMSEY.R <- array(NA, dim=c(nrows,ncols))
for (i in 1:ncols) {
    this.Lmeas <- x.RAMSEY.all[[i]]$Lmeas
    all.RAMSEY.Lmeas[1:length(this.Lmeas),i] <- this.Lmeas
    this.R <- x.RAMSEY.all[[i]]$R
    all.RAMSEY.R[1:length(this.R),i] <- this.R
}
mean.RAMSEY.Lmeas <- apply(all.RAMSEY.Lmeas, 1, function(x) mean(x, na.rm=T))
mean.RAMSEY.R <- apply(all.RAMSEY.R, 1, function(x) mean(x, na.rm=T))

par(mfrow=c(1,2), oma=rep(0.5,4), mar=c(5,5,1,1))
plot.new()
plot.window(xlim=c(0,50), ylim=c(0,3))
axis(1)
axis(2)
box('plot')
mtext(side=1, line=3, 'Age')
mtext(side=2, line=3, 'Length')
points(mean.2c.age, mean.2c.Lmeas, type='l', lwd=1.5, col='blue')
points(mean.4b.age, mean.4b.Lmeas, type='l', lwd=1.5, col='green')
points(mean.GLEN.age, mean.GLEN.Lmeas, type='l', lwd=1.5, col='gold')
points(mean.RAMSEY.age, mean.RAMSEY.Lmeas, type='l', lwd=1.5, col='red')

plot.new()
plot.window(xlim=c(0,50), ylim=c(0,30))
axis(1)
axis(2)
box('plot')
mtext(side=1, line=3, 'Age')
mtext(side=2, line=3, 'Cum. Eggs')
points(mean.2c.age, mean.2c.R, type='l', lwd=1.5, col='blue')
points(mean.4b.age, mean.4b.R, type='l', lwd=1.5, col='green')
points(mean.GLEN.age, mean.GLEN.R, type='l', lwd=1.5, col='gold')
points(mean.RAMSEY.age, mean.RAMSEY.R, type='l', lwd=1.5, col='red')

dev.copy2pdf(file='Bill_mean_fit_trajectories.pdf')

## Calculate the median trajectories for each genotype
## Determine the length of the longest-lived individual
nrows <- max(sapply(x.2c.all, function(x) nrow(x)))
ncols <- length(x.2c.all)
median.2c.age <- seq(0,(nrows-1)*0.1,0.1)
all.2c.Lmeas <- array(NA, dim=c(nrows,ncols))
all.2c.R <- array(NA, dim=c(nrows,ncols))
for (i in 1:ncols) {
    this.Lmeas <- x.2c.all[[i]]$Lmeas
    all.2c.Lmeas[1:length(this.Lmeas),i] <- this.Lmeas
    this.R <- x.2c.all[[i]]$R
    all.2c.R[1:length(this.R),i] <- this.R
}
median.2c.Lmeas <- apply(all.2c.Lmeas, 1, function(x) median(x, na.rm=T))
median.2c.R <- apply(all.2c.R, 1, function(x) median(x, na.rm=T))

nrows <- max(sapply(x.4b.all, function(x) nrow(x)))
ncols <- length(x.4b.all)
median.4b.age <- seq(0,(nrows-1)*0.1,0.1)
all.4b.Lmeas <- array(NA, dim=c(nrows,ncols))
all.4b.R <- array(NA, dim=c(nrows,ncols))
for (i in 1:ncols) {
    this.Lmeas <- x.4b.all[[i]]$Lmeas
    all.4b.Lmeas[1:length(this.Lmeas),i] <- this.Lmeas
    this.R <- x.4b.all[[i]]$R
    all.4b.R[1:length(this.R),i] <- this.R
}
median.4b.Lmeas <- apply(all.4b.Lmeas, 1, function(x) median(x, na.rm=T))
median.4b.R <- apply(all.4b.R, 1, function(x) median(x, na.rm=T))

nrows <- max(sapply(x.GLEN.all, function(x) nrow(x)))
ncols <- length(x.GLEN.all)
median.GLEN.age <- seq(0,(nrows-1)*0.1,0.1)
all.GLEN.Lmeas <- array(NA, dim=c(nrows,ncols))
all.GLEN.R <- array(NA, dim=c(nrows,ncols))
for (i in 1:ncols) {
    this.Lmeas <- x.GLEN.all[[i]]$Lmeas
    all.GLEN.Lmeas[1:length(this.Lmeas),i] <- this.Lmeas
    this.R <- x.GLEN.all[[i]]$R
    all.GLEN.R[1:length(this.R),i] <- this.R
}
median.GLEN.Lmeas <- apply(all.GLEN.Lmeas, 1, function(x) median(x, na.rm=T))
median.GLEN.R <- apply(all.GLEN.R, 1, function(x) median(x, na.rm=T))

nrows <- max(sapply(x.RAMSEY.all, function(x) nrow(x)))
ncols <- length(x.RAMSEY.all)
median.RAMSEY.age <- seq(0,(nrows-1)*0.1,0.1)
all.RAMSEY.Lmeas <- array(NA, dim=c(nrows,ncols))
all.RAMSEY.R <- array(NA, dim=c(nrows,ncols))
for (i in 1:ncols) {
    this.Lmeas <- x.RAMSEY.all[[i]]$Lmeas
    all.RAMSEY.Lmeas[1:length(this.Lmeas),i] <- this.Lmeas
    this.R <- x.RAMSEY.all[[i]]$R
    all.RAMSEY.R[1:length(this.R),i] <- this.R
}
median.RAMSEY.Lmeas <- apply(all.RAMSEY.Lmeas, 1, function(x) median(x, na.rm=T))
median.RAMSEY.R <- apply(all.RAMSEY.R, 1, function(x) median(x, na.rm=T))

par(mfrow=c(1,2), oma=rep(0.5,4), mar=c(5,5,1,1))
plot.new()
plot.window(xlim=c(0,50), ylim=c(0,3))
axis(1)
axis(2)
box('plot')
mtext(side=1, line=3, 'Age')
mtext(side=2, line=3, 'Length')
points(median.2c.age, median.2c.Lmeas, type='l', lwd=1.5, col='blue')
points(median.4b.age, median.4b.Lmeas, type='l', lwd=1.5, col='green')
points(median.GLEN.age, median.GLEN.Lmeas, type='l', lwd=1.5, col='gold')
points(median.RAMSEY.age, median.RAMSEY.Lmeas, type='l', lwd=1.5, col='red')

plot.new()
plot.window(xlim=c(0,50), ylim=c(0,30))
axis(1)
axis(2)
box('plot')
mtext(side=1, line=3, 'Age')
mtext(side=2, line=3, 'Cum. Eggs')
points(median.2c.age, median.2c.R, type='l', lwd=1.5, col='blue')
points(median.4b.age, median.4b.R, type='l', lwd=1.5, col='green')
points(median.GLEN.age, median.GLEN.R, type='l', lwd=1.5, col='gold')
points(median.RAMSEY.age, median.RAMSEY.R, type='l', lwd=1.5, col='red')

dev.copy2pdf(file='Bill_median_fit_trajectories.pdf')

## Try to recreate something that can better compare against Fig. 2 in
## Adriana and Bill's paper. Specifically, make one plot of all of the
## asymptotic lengths and egg rates.
all.2c.Lmax <- apply(all.2c.Lmeas, 2, function(x) max(x, na.rm=T))
all.4b.Lmax <- apply(all.4b.Lmeas, 2, function(x) max(x, na.rm=T))
all.GLEN.Lmax <- apply(all.GLEN.Lmeas, 2, function(x) max(x, na.rm=T))
all.RAMSEY.Lmax <- apply(all.RAMSEY.Lmeas, 2, function(x) max(x, na.rm=T))

## Turn cumulative eggs into a calculation of mean daily egg rate for each individual
all.2c.eggrate <- apply(all.2c.R, 2, function(x) max(x, na.rm=T))/pars.2c$aod
all.4b.eggrate <- apply(all.4b.R, 2, function(x) max(x, na.rm=T))/pars.4b$aod
all.GLEN.eggrate <- apply(all.GLEN.R, 2, function(x) max(x, na.rm=T))/pars.GLEN$aod
all.RAMSEY.eggrate <- apply(all.RAMSEY.R, 2, function(x) max(x, na.rm=T))/pars.RAMSEY$aod

par(mfrow=c(1,2), oma=rep(0.5,4), mar=c(5,5,1,1))
plot.new()
plot.window(xlim=c(0,1), ylim=c(2,3))
axis(2)
box('plot')
mtext(side=2, line=3, 'Lmax')
points(rep(0.2,length(all.2c.Lmax)), all.2c.Lmax, col='blue')
points(rep(0.4,length(all.4b.Lmax)), all.4b.Lmax, col='green')
points(rep(0.6,length(all.GLEN.Lmax)), all.GLEN.Lmax, col='gold')
points(rep(0.8,length(all.RAMSEY.Lmax)), all.RAMSEY.Lmax, col='red')
abline(h=median(all.2c.Lmax), col='blue')
abline(h=median(all.4b.Lmax), col='green')
abline(h=median(all.GLEN.Lmax), col='gold')
abline(h=median(all.RAMSEY.Lmax), col='red')

plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1))
axis(2)
box('plot')
mtext(side=2, line=3, 'Mean egg rate')
points(rep(0.2,length(all.2c.eggrate)), all.2c.eggrate, col='blue')
points(rep(0.4,length(all.4b.eggrate)), all.4b.eggrate, col='green')
points(rep(0.6,length(all.GLEN.eggrate)), all.GLEN.eggrate, col='gold')
points(rep(0.8,length(all.RAMSEY.eggrate)), all.RAMSEY.eggrate, col='red')
abline(h=median(all.2c.eggrate), col='blue')
abline(h=median(all.4b.eggrate), col='green')
abline(h=median(all.GLEN.eggrate), col='gold')
abline(h=median(all.RAMSEY.eggrate), col='red')

dev.copy2pdf(file='Bill_Lmax_eggrate_plot.pdf')
