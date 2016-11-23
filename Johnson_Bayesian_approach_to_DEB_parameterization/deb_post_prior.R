## In this file we put the functions that calculate prior and
## posterior probabilities of the parameters. The function that
## returns the posterior probabilities must be called
## "log.post.params", and the prior functions should in
## "log.prior.params" (which will only be called from
## log.post.params). The user can update these functions so they are
## appropriate for their data and parameters. Additionally the
## function "make.hypers" makes an appropriately structured set of
## hyper parameters (which determine the prior) for the deb.mcmc
## function.



## the following functions define the posterior and prior
## distributions for the parameters of interest: J.EAm, y.EX, v, k.M,
## g, kap, k.J, M.HP, shape, gamma, X.h, vol. There is also a function
## to automatically generate the hyper parameters for the priors

log.post.params<-function(samps, data, samp.p, sim.data){

  e.c<-0.0001

  ## observation model
  L.M<-as.numeric(samps["L.M"])
  gamma<-as.numeric(samps["gamma"])
  
  l.temp<-sim.data$l*L.M
  n.temp<-sim.data$n*gamma

  sd.L<-as.numeric(samps["sd.L"])
  sd.Negg<-as.numeric(samps["sd.Negg"])


  llik.L<-sum(dlnorm(data$L, meanlog=log(l.temp+e.c), sdlog=sd.L, log=TRUE))
  llik.Negg<-sum(dlnorm(data$Negg+0.000001, meanlog=log(n.temp+e.c), sdlog=sd.Negg, 
log=TRUE))
   
  llik<-llik.L+llik.Negg
  
  lprior<-sum(log.prior.params(samps, samp.p))

  if(is.na(llik)|!is.finite(llik)) stop("something is wrong in the likelihood bit")
  if(is.na(lprior)|!is.finite(lprior)) stop("something is wrong in the prior bit")

  return( llik + lprior )	
}



## eventually change this to make it more efficient
log.prior.params<-function(samps, samp.p){
  
  lp<-NULL
  cnt<-0
  len<-length(samp.p)
  
  ##print(c(as.numeric(samp), w.p, hyper[1]))
  ## this won't work as it is for joint proposals
  for(i in 1:len){
    l<-length(samp.p[[i]]$params)
    for(j in 1:l){
      cnt<-cnt+1
      p<-samp.p[[i]]$params[j]
      s<-as.numeric(samps[p])

      if(l==1){
        param1<-samp.p[[i]]$hyp[1]
        param2<-samp.p[[i]]$hyp[2]
      }
      else{
        param1<-samp.p[[i]]$hyp[j,1]
        param2<-samp.p[[i]]$hyp[j,2]
      }
      
      if( p == "J.EAm" ) lp<-c(lp, J.EAm=dgamma(s, shape=param1, rate=param2, log=TRUE))
      if( p == "y.EX" ) lp<-c(lp, y.EX=dgamma(s, shape=param1, rate=param2, log=TRUE))
      if( p == "v" ) lp<-c(lp, v=dgamma(s, shape=param1, rate=param2, log=TRUE))
      if( p == "k.M" ) lp<-c(lp, k.M=dgamma(s, shape=param1, rate=param2, log=TRUE))
      if( p == "g" ) lp<-c(lp, g=dgamma(s, shape=param1, rate=param2, log=TRUE))
      if( p == "kap" ) lp<-c(lp, kap=dbeta(s, shape1=param1, shape2=param2,  log=TRUE))
      if( p == "k.J" ) lp<-c(lp, k.J=dgamma(s, shape=param1, rate=param2, log=TRUE))
      if( p == "M.HP" ) lp<-c(lp, M.HP=dgamma(s, shape=param1, rate=param2,  log=TRUE))
      if( p == "shape" ) lp<-c(lp, shape=dgamma(s, shape=param1, rate=param2, log=TRUE))
      if( p == "gamma" ) lp<-c(lp, gamma=dnorm(s, mean=param1, sd=param2, log=TRUE))
      if( p == "X.h" )  lp<-c(lp, X.h=dgamma(s, shape=param1, rate=param2, log=TRUE))
      if( p == "vol" )  lp<-c(lp, vol=dgamma(s, shape=param1, rate=param2, log=TRUE))
      if( p == "sd.L" )  lp<-c(lp, sd.L=dgamma(s, shape=param1, rate=param2, log=TRUE))
      ##if( p == "sd.L" )  lp<-c(lp, sd.L=dunif(s, param1, param2, log=TRUE))
      if( p == "sd.Negg" )  lp<-c(lp, sd.Negg=dgamma(s, shape=param1, rate=param2, log=TRUE))
      ##if( p == "sd.Negg" )  lp<-c(lp, sd.Negg=dunif(s, param1, param2, log=TRUE))
      if( p == "l.M.HP" ) lp<-c(lp, l.M.HP=dnorm(s, mean=param1, sd=param2,  log=TRUE))
      if( p == "L.M" ) lp<-c(lp, L.M=dgamma(s, shape=param1, rate=param2,  log=TRUE))
      if( p == "F.eff" ) lp<-c(lp, F.eff=dgamma(s, shape=param1, rate=param2,  log=TRUE))

     if(!is.finite(lp[cnt])) stop("one of the prior probs isn't finite") 
    }
  }

  return(lp)

}


## here's a function to make the hyper params in the correct format

make.hypers<-function(J.EAm=c(1,1), y.EX=c(1,1), v=c(1,1), k.M=c(1,1),
                      g=c(1,1), kap=c(2,2),  k.J=c(1,1), M.HP=c(1,5),
                      shape=c(1,1), gamma=c(500,100), X.h=c(1,1),
                      vol=c(1,1), sd.L=c(1,1), sd.Negg=c(1,1), l.M.HP=c(-1,5),
                      L.M=c(1,1), F.eff=c(1,1)){

  hyper<-list(J.EAm=J.EAm, y.EX=y.EX, v=v, k.M=k.M, g=g, kap=kap,
              k.J=k.J, M.HP=M.HP, shape=shape, gamma=gamma,
              X.h=X.h, vol=vol, sd.L=sd.L, sd.Negg=sd.Negg, l.M.HP=l.M.HP,
              L.M=L.M, F.eff=F.eff)
  
  return(hyper)

}


