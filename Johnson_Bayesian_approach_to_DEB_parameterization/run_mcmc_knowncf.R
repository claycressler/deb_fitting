## here's code to run a sample of our nifty mcmc
## run with: nohup R CMD batch run_mcmc.R &



source("deb_mcmc.R")
source("deb_model.R")
source("deb_solver.R")

## need to source this one after deb_mcmc.R to tell it to use this
## versio of the prior/posterior calcs instead of the versions in
## deb_post_prior.R
source("deb_post_prior.R")

##library("coda")

set.seed(12345)
plt=TRUE
cnt<-10
N<-300
my.par<-c(4,2)


w<-c("sd.L", "sd.Negg", "l.M.HP", "L.M", "v", "kap", "g") ##, "F.eff")
start<-c(sd.L=0.1, sd.Negg=0.1, F.eff=0.35, g=0.498, kap=0.1, v=1.48,
         l.M.HP=log(0.00575), L.M=4)
samp.type<-list(sd.L="rw", sd.Negg="rw", F.eff="ind", g="rw", kap="rw",
                v="rw", l.M.HP="rw", L.M="rw", X.h="ind")
prop.var<-list(sd.L=c(2,3), sd.Negg=c(2,3), F.eff=0.0001, g=0.006, L.M=0.006,
               kap=0.004, v=0.017, l.M.HP=0.001)
hyper<-NULL


inits<-inits.sim<-setinits.DEB(X0=0.2)
params.sim<-params<-setparams.DEB(var.food=0, feed.space=NULL, feed.level=NULL)

inits[2]<-inits.sim[2]<-as.numeric(PHI(params.sim["X.h"], inits.sim[1]))
sw.name<-NULL
mp.name<-NULL
Tmax<-50


ss<-0.01
data<-solve.DEB(DEB.daphnia2, params.sim, inits.sim, Tmax=Tmax, numsteps=NULL,
                which="cont", sizestep=ss, myswitch=sw.name, mymap=mp.name)

if(plt) plot.DEB(data, scaled.length=TRUE, scale=10)
data<-make.obs(data, params, w.t=1, Tmax=Tmax)

if(plt) plot.DEB.red(data)

set.seed(1)

mcmc.p<-setup.mcmc.params(w.p=w, samp.type=samp.type, start=start,
                          prop.var=prop.var, joint=FALSE, n.joint=1, scale.var=7)

samps<-deb.mcmc(N, data, params, inits, sim=DEB.daphnia2, mcmc.p, Tmax,
                cnt=cnt, burnin=0, plot=plt, sizestep=ss, which="cont",
                my.par=my.par,  myswitch=NULL, mymap=NULL)

samps<-samps$samps
if(is.null(hyper)) hyper<-make.hypers()
w<-c("sd.L", "sd.Negg", "v", "l.M.HP", "kap", "g", "L.M")##, "F.eff")


save(data, mcmc.p, params, inits, hyper, samps, w, file="test_known_cf.Rsave")

