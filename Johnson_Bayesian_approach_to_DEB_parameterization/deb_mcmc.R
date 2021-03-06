

library("mvtnorm")

## functions for determining the posterior and prior probs are
## contained here:
## source("deb_post_prior.R")


## Bayesian inference for a deterministic DEB model (with models
## solved via an ode solver in the function specified with argument
## "sim") with an observation model, and with observation error
## variances specified in sds.

deb.mcmc<-function(N, data, all.params, inits, sim=DEB1, samp.p,
                   Tmax, cnt=10, burnin=0.1,
                   plot=TRUE, sizestep=0.01, w.t=1, which=2,
                   test=TRUE, my.par=c(1,1),  myswitch=NULL,
                   mymap=NULL)
{

  ## first we calculate a few lengths, etc, that we use for the for
  ## loops later. 
  
  l<-length(samp.p)
  np<-length(all.params)
  p<-NULL

  ltot<-0
  for(i in 1:l) ltot<-ltot+length(samp.p[[i]]$params)
  
  ## for testing, here is code that calculates (and prints out) the
  ## posterior prob of the "real" parameters, which can be passed in
  ## through all.params
  
  if(test){
    sim.old<-make.states(sim, all.params, inits, Tmax, which=which, sizestep, w.t,
                          myswitch=myswitch, mymap=mymap)
    prob.old<-log.post.params(all.params, data, samp.p, sim.old)
    print(paste("posterior probability of the real parameters= ", prob.old, sep=""))
  }

  ## build a data frame to hold the posterior samples. These include
  ## "samples" even for the parameters that are being held fixed, to
  ## make the code more straightforward and generic, as well as being
  ## useful for debugging. The samps structure will also keep track of 
  ## the log posterior probability of that particular sample.

  samps<-data.frame(matrix(0, ncol=np+1, nrow=N))
  names(samps)<-c(names(all.params), "lpost")

  ## initialize the samps structure and p
  samps[1,1:np]<-p<-all.params  
  
  ## initialize the data frame with the starting values passed into
  ## the function in the structure samp.p
  for(i in 1:l){
    samps[1, samp.p[[i]]$params]<-samp.p[[i]]$start
    p[samp.p[[i]]$params]<-samp.p[[i]]$start
  }


  ## run the data simulation to make the underlying states (for
  ## determining the likelihoods) using the parameters stored in p.
  
  sim.old<-make.states(sim, p, inits, Tmax, which=which, sizestep, w.t,
                        myswitch=myswitch, mymap=mymap)

  ## check the posterior probability to make sure you have reasonable
  ## starting values, and to initialize prob.old
  samps$lpost[1]<-log.post.params(samps[1,], data, samp.p, sim.old)
  
  print(paste("initial posterior probability = ", samps$lpost[1], sep=""))
  
  if(!is.finite(samps$lpost[1])) stop("bad starting values")

  ## now we begin the MCMC

  out<-list(s=samps[1,], p=p, sim.old=sim.old)
  
  for(i in 2:N){
    ## printing and plotting output so we can watch the progress
    if(i%%cnt == 0){
      print(paste("sample number", i, sep=" "))
      if(plot) plot.output(i, samps, samp.p, l, ltot, my.par)
    }

    ## the meat of the MCMC is found in the function update.samps (see below)

    out<-update.sample(samps[i-1,], samp.p, data, sim, inits, out,
                       Tmax, sizestep, w.t, l, which, i, cnt,
                        myswitch=myswitch, mymap=mymap)
    samps[i,]<-out$s
    if(test){
      if(-samps$lpost[i-1]+samps$lpost[i-1]<=-10){
        stop("we've had a really large swing in the posterior prob")
      }
    }
    
  }
	
  ##plot(b[100:N],type="l")
  lim<-min(1, round(burnin*N))
  samps <- samps[lim:N,]

  return(list(samps=samps))
	
}



update.sample<-function(samps, samp.p, data, sim, inits, out, Tmax, sizestep,
                        w.t, l, which, i, cnt,  myswitch=NULL, mymap=NULL, test=TRUE)
{
  ## read in some bits
  s<-samps
  sim.old<-out$sim.old
  p<-out$p

  x<-1:l
  s.x<-sample(x)
  
  for(k in s.x){   

    s.new<-s
    p.new<-p
   
    q<-propose.params(s, samp.p[[k]])

    ## automatically reject if the proposed parameter value is
    ## outside of the reasonable limits (i.e. < 0 )
    zeros<-0
    zeros<-check.zeros(samp.p[[k]], q$b)
    if(zeros){
      print("proposed a param < 0. moving on")
      next
    }
    ## write the proposed params into the p.new and s.new.
    
    for(j in 1:length(samp.p[[k]]$params)){
      ww<-samp.p[[k]]$params[j]
      p.new[ww]<-s.new[ww]<-q$b[j]
    }

    ## simulate the dynamics forward with the new parameters
    sim.new<-make.states(sim, p.new, inits, Tmax, which=which, sizestep, w.t,
                          myswitch=myswitch, mymap=mymap)
 
    ## The posteriorprob of the previous sample is saved as
    ## s$lpost. If we accept a draw, we will set s$lpost<-s.new$lpost
    
    s.new$lpost<-log.post.params(s.new, data, samp.p, sim.new)
    
    if(is.finite(s.new$lpost) && is.finite(s$lpost)){
      A<-exp( s.new$lpost + q$lbak - s$lpost - q$lfwd )
    }
    else{
      A<-0
      print("whoops! must have proposed outside the correct region")
    }

    ## print some output so we can follow the progress
    if(i%%cnt==0){
      print(paste("proposing " , samp.p[[k]]$params, ": prob.old = ",
                  signif(s$lpost, digits=5),
                  "; prob.new = ", signif(s.new$lpost, digits=5), "; A = ",
                  signif(A, digits=5),
                  sep=""))
    }
    
    ## take a draw from a unif distrib, and use it to accept/reject
    u<-runif(1)    
    if( u < A ){ ## accept
      sim.old<-sim.new
      p<-p.new
      s<-s.new
    }
  }

  return(list(s=s, p=p, sim.old=sim.old))
  
}


check.zeros<-function(s.p, q.b){
  z<-0
  for(j in 1:length(s.p$params)){
    if(s.p$params[j]=="l.M.HP") next
    if(q.b[j]<0)z<-1
    if(s.p$params[j]=="kap"){
      if(q.b[j]>1) z<-1
    }
  }
  return(z)
}




propose.params<-function(samps, s.p)
{
  if(length(s.p$params)==1){
    ##print(paste(s.p$params, " proposing single ", sep=" "))
    q<-propose.single(samps, s.p)
  }
  else{
    ##print(paste(s.p$params, " proposing jointly ", sep=" "))
    q<-propose.joint(samps, s.p)
  }         
  return(q)
     
}

## I'm feeding in the variance, so I need to take the square root....
propose.single<-function(samps, s.p)##, i, freq=50, size=50 )##l=5, h=6)
{

  b<-as.numeric(samps[s.p$params])
  var<-s.p$var
  type<-s.p$type
  hyps<-s.p$hyp
  
  if(type=="rw"){
    if(length(var)>1){
      l<-var[1]
      h<-var[2]
      b.new<-runif(1, l/h*b, h/l*b)
      lfwd<-dunif(b.new, l/h*b, h/l*b, log=TRUE)
      lbak<-dunif(b, l/h*b.new, h/l*b.new, log=TRUE)
    }
    else{
      sd<-sqrt(var)
      b.new<-rnorm(1, b, sd=sd)
      lfwd<-dnorm(b.new, b, sd=sd, log=TRUE)
      lbak<-dnorm(b, b.new, sd=sd, log=TRUE)
      
    }
    return(list(b=b.new, lfwd=lfwd, lbak=lbak))
  }
  else if(type=="ind"){
    out<-prior.draw(b, hyps, s.p$params)
    return(out)
  }
  
}	



propose.joint<-function(samp, samp.p){

  
  b<-NULL
  if(samp.p$type=="rw"){
    b<-as.numeric(samp[samp.p$params])
    sigma<-samp.p$var
    
    b.new<-rmvnorm(1, mean=b, sigma=sigma)
    lfwd<-dmvnorm(b.new, b, sigma, log=TRUE)
    lbak<-dmvnorm(b, b.new, sigma, log=TRUE)
  }
  else if(samp.p$type=="ind"){
    if(is.null(samp.p$mean)) stop("not enough info for the independence sampler")
    mean<-as.numeric(samp.p$mean)
    b<-as.numeric(samp[samp.p$params])
    sigma<-samp.p$var
    
    b.new<-rmvnorm(1, mean=mean, sigma=sigma)
    lfwd<-dmvnorm(b.new, mean, sigma, log=TRUE)
    lbak<-dmvnorm(b, mean, sigma, log=TRUE)
  }
  
  ##print(c(b, b.new, lfwd, lbak))


  ##samp[s]<-b.new

  ##stop()
  return(list(b=b.new, lfwd=lfwd, lbak=lbak))

}



plot.output<-function(i, samps, samp.p, l, ltot, my.par=c(2,4), plot.post=TRUE){
  
  if( ltot > 1 ) par(mfrow=my.par, bty="n")
  for( j in 1:l ){
    ww<-samp.p[[j]]$params
    for(k in 1:length(ww)){
      plot(samps[1:(i-1),ww[k]], type="l", xlab="sample", ylab=ww[k])
    }
  }
  if(plot.post){
    my.min<-20
    if(i>my.min){
      x<-seq(my.min, (i-1), by=1)
      plot(x, samps$lpost[x], type="l", xlab="sample", ylab="log posterior prob")
    }
  }
}


prior.draw<-function(b, hyp, p){
  
  param1<-hyp[1]
  param2<-hyp[2]
 
  gams<-c("J.EAm", "y,EX", "v", "k.M", "g", "k.J", "M.HP", "shape", "X.h",
          "vol", "L.M", "F.eff") ##, "sd.L", "sd.Negg"
  norms<-c("l.M.HP", "gamma")
  betas<-c("kap")
  unifs<-c("sd.L", "sd.Negg")
      
  if( p %in% gams ){
    b.new<-rgamma(1, shape=param1, rate=param2)
    lfwd<-dgamma(b.new, shape=param1, rate=param2, log=TRUE)
    lbak<-dgamma(b, shape=param1, rate=param2, log=TRUE)
  }
  else if( p %in% norms ){
    b.new<-rnorm(1, mean=param1, sd=param2)
    lfwd<-dnorm(b.new, mean=param1, sd=param2, log=TRUE)
    lbak<-dnorm(b, mean=param1, sd=param2, log=TRUE)
  }
  else if( p %in% betas ){
    b.new<-rbeta(1,  shape1=param1, shape2=param2)
    lfwd<-dbeta(b.new,  shape1=param1, shape2=param2, log=TRUE)
    lbak<-dbeta(b, shape1=param1, shape2=param2, log=TRUE)
  }
  else if( p %in% unifs ){
    b.new<-runif(1,  param1, param2)
    lfwd<-dunif(b.new,  param1, param2, log=TRUE)
    lbak<-dunif(b, param1, param2, log=TRUE)
  }

  return(list(b=b.new, lfwd=lfwd, lbak=lbak)) 
}


