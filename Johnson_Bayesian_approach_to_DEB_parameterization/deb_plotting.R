source("deb_mcmc.R")
source("deb_solver.R")
source("deb_model.R")

## need to source this one after deb_mcmc.R to tell it to use this
## version of the prior/posterior calcs instead of the versions in
## deb_post_prior.R
source("deb_post_prior.R")


library(IDPmisc)
  

plot.samps<-function(samps, w, burnin, N, col=TRUE, thin=1){
  thinned<-seq(burnin, N, by=thin)
  if(col) ipairs(samps[thinned,w], ztransf = function(x){x[x<1] <- 1; log2(x)})
  else     plot(samps[thinned,w])
}


plot.trace<-function(samps, w, burnin, N, my.par, LL=TRUE, thin=1){
  
  thinned<-seq(burnin, N, by=thin)
  
  par(mfrow=my.par, bty="n")
  for(i in w) plot(thinned, samps[thinned,i], type="l", xlab="sample", ylab=i)
  if(LL) plot(thinned, samps$lpost[thinned], type="l", xlab="sample", ylab="log posterior")
  
}


plot.marg<-function(samps, w, hyper, burnin, N, my.par, thin=1, true.params=NULL){
  thinned<-seq(burnin, N, by=thin)
  
  par(mfrow=my.par, bty="n")
  x<-seq(0,100, length=5000)
  for(i in w){
    h<-unlist(hyper[i])
    hist(samps[thinned,i], freq=FALSE, main=" ", xlab=i)
    if(i=="kap"){
      lines(x, dbeta(x, shape1=as.numeric(h[1]),
                     shape2=as.numeric(h[2])), col=2)
    }
    else if(i=="l.M.HP"){
      y<-seq(-20,20,length=5000)
      lines(y, dnorm(y, mean=as.numeric(h[1]),
                     sd=as.numeric(h[2])), col=2)
    }
    else{
      lines(x, dgamma(x, shape=as.numeric(h[1]), rate=as.numeric(h[2])), col=2)
    }
    points(mean(samps[thinned,i]), 0, col="red", pch=17, cex=2.5)
    if(!is.null(true.params)) points(true.params[i], 0, col="blue", pch=19, cex=2)    
  }
  hist(samps$v[thinned]/(samps$g[thinned]*samps$L.M[thinned])*0.1, main="", xlab="k.J",
       freq=FALSE, breaks=20)
  if(!is.null(true.params)) points(true.params["v"]/(true.params["g"]*true.params["L.M"])*0.1, 0,
                                   pch=19, col="blue", cex=2)
  points(mean(samps$v[thinned]/(samps$g[thinned]*samps$L.M[thinned])*0.1),0, pch=17,
         col="red", cex=2.5)

}


deb.make.sims<-function(sim=DEB.daphnia2, samps, inits, thinned, Tmax, ss,
                    which="cont", myswitch=NULL, mymap=NULL){

  out<-list()
  data.sim<-NULL
  for(i in 1:length(thinned)){
    prms<-mean(samps[thinned[i],])
    data<-solve.DEB(sim, prms, inits, Tmax=Tmax, numsteps=NULL,
                which=which, sizestep=ss, myswitch=myswitch, mymap=mymap,
                    adj.length=TRUE)
    names(data)<-c("t", "F", "e", "l", "M", "R")
    out[[i]]<-data
  }
  return(out)
}

deb.make.states<-function(sim.data, samps, thinned, ss){
  ll<-length(sim.data[[1]]$l)
  L<-matrix(NA, nrow=length(thinned), ncol=ll)
  Negg<-matrix(NA, nrow=length(thinned), ncol=ll)
  for(i in 1:length(thinned)){
    L[i,]<-sim.data[[i]]$l*samps$L.M[thinned[i]]
    Negg[i,]<-sim.data[[i]]$R*samps$gamma[thinned[i]]
  }
  return(list(L=L, Negg=Negg))
}


deb.sim.quants<-function(sim.data, l, byCol=TRUE,
                     probs=c(0.025, 0.975)){

  q<-matrix(NA, nrow=length(probs), ncol=l)
  if(byCol) for(i in 1:l) q[,i]<-quantile(sim.data[,i], probs, na.rm=TRUE)
  else for(i in 1:l) q[i,]<-quantile(sim.data[i,], probs, na.rm=TRUE)
  
  return(q)
}


deb.sim.quants.LN<-function(sim.data, samps.sd, l, N.sims, byCol=TRUE, probs=c(0.025, 0.975)){

  ec<-0.0001

  q.low<-q.high<-matrix(NA, nrow=N.sims, ncol=l)
  q<-matrix(NA, nrow=length(probs), ncol=l)

  q.temp<-0
  for(i in 1:l){
    q.low[,i]<-qlnorm(probs[1], meanlog=log(sim.data[,i]+ec), sdlog=samps.sd)
    q.high[,i]<-qlnorm(probs[2], meanlog=log(sim.data[,i]+ec), sdlog=samps.sd)
  }

  if(byCol){
    for(i in 1:l){
      q[1,i]<-quantile(q.low[,i], probs[1], na.rm=TRUE)
      q[2,i]<-quantile(q.high[,i], probs[2], na.rm=TRUE)
    }
  }
  else{
    for(i in 1:l){
      q[i,1]<-quantile(q.low[i,], probs[1], na.rm=TRUE)
      q[i,2]<-quantile(q.high[i,], probs[2], na.rm=TRUE)
    }
  }
  ##stop()
  return(q)
  

}

add.sim.lines<-function(t, sim.data=NULL, q=NULL, q2=NULL, mycol="blue", lwd=1){

  if(!is.null(sim.data)){
    ## plot the predicted mean
    if(!is.null(dim(sim.data)[1])) lines(t, colMeans(sim.data, na.rm=TRUE), col=mycol, lwd=lwd)
    else lines(t, sim.data, col=mycol, lwd=lwd)
    
    if(!is.null(q)){
      ## plot the predicted quantiles for the mean
      lines(t, q[1,], col=mycol, lty="dashed", lwd=(lwd+1))
      lines(t, q[2,], col=mycol, lty="dashed", lwd=(lwd+1))
    }
    if(!is.null(q2)){
      ## plot the predicted quantiles for the mean
      lines(t, q2[1,], col=mycol, lty="dotted", lwd=(lwd+2))
      lines(t, q2[2,], col=mycol, lty="dotted", lwd=(lwd+2))
    }
    
  }
}

makeplot2<-function(data, data.obs, sim.data=NULL, q=NULL, q2=NULL, myxlim=c(0,8),
                    mycol="blue", add=FALSE){

  if(!add){
    boxplot((data$counts+1)~data$time, ylab="Zoospores",
         ylim=c(1, 1000), xlim=myxlim, log="y",
         yaxt = "n", xlab = "Time (d)", col=mycol)
       axis(2, c(1,10,100,1000), labels= c(0,10,100,1000), las=1)
  }
  else boxplot((data$counts+1)~data$time, col=mycol, add=TRUE)
  for(i in 1:n.obs) points(data.obs$t, (data.obs$Z[i,]+1), col=mycol)

  if(!is.null(sim.data)){
    ## plot the predicted mean
    lines(sim.data$t, colMeans(sim.data$Z)+1, col=mycol)
    
    if(!is.null(q)){
      ## plot the predicted quantiles for the mean
      lines(sim.data$t, q[1,]+1, col=mycol, lty="dashed")
      lines(sim.data$t, q[2,]+1, col=mycol, lty="dashed")
    }
    if(!is.null(q2)){
      ## plot the predicted quantiles for the mean
      lines(sim.data$t, q2[1,]+1, col=mycol, lty="dotted")
      lines(sim.data$t, q2[2,]+1, col=mycol, lty="dotted")
    }
  }
}




