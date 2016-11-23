##install PBSddesolve first
library(Matrix)
library(PBSddesolve)

## Daphnia model (see Nisbet et al 2010 PTRS).

DEB.daphnia2<-function(t,y,p){

   J.EAm <- p["J.EAm"]
   y.EX <- p["y.EX"]
   v <- p["v"]
   L.M <- p["L.M"]
   g <- p["g"]
   kap <- p["kap"]
   k.J <- p["k.J.scale"]*(v/(g*L.M))
   M.HP <- exp(p["l.M.HP"])
   gamma <- p["gamma"]
   X.h <- p["X.h"]
   vol <- p["vol"]
   fb.f<-p["var.food"] ## this now determines whether or not there is
                      ## feedback in the food, not if the food is
                      ## constant 

   x<-y[1]
   e<-y[2]
   l<-y[3]
   M.H<-y[4]
   M.R<-y[5]
   
   phi<-PHI(X.h, x) #y[1]/(y[1]+X.h)
   
   F<-food(t,X.h,x)   
   
   ## here Delta = -(Delta in the paper) -- ie, if >0 sufficient reserve
   ## to pay somatic maintenance
   Delta <- - L.M^2*J.EAm* kap * g * (l-e)*l^2/(e + g)
   ##stop(paste("delta=", Delta, sep=""))

   ## p.c, or nearly  
   p.c<-(1 - kap) * l^2 * L.M^2 * e * J.EAm /(e+g)
   
   if(Delta > 0){
     dy1 <- -J.EAm*L.M^2/ (vol * y.EX) * l^2 * phi * fb.f + F #dX
     dy2 <- v / L.M * (phi-e)/l #de
     dy3 <- v/(3*L.M) * (e-l)/(e + g)#dl
     
     ## Development and reproduction
     if(y[4] < M.HP){
       dy4 <-  p.c  - k.J * M.H #dMH
       dy5 <- 0 # dMR
     }
     else{
       dy4 <- 0
       dy5 <-  p.c - k.J * M.H
     }
   }
   else{## Delta>0, i.e. J.EC > J.EJ + J.ES
     dy1 <- -J.EAm*L.M^2/ (vol * y.EX) * l^2 * phi * fb.f + F #dX
     dy2 <- v / L.M * (phi-e)/l #de
     dy3 <- 0 #dl
     ## Development and reproduction
     if(y[4] < M.HP){
       dy4 <- p.c - k.J * M.H + Delta #dMH
       dy5 <- 0  # dMR
     }
     else{
       dy4 <- 0
       dy5 <- p.c - k.J * M.H + Delta
     }
     ##}
     if(dy4<0 | dy5<0){# death
       dy1 <- F
       dy2<- dy3 <- dy4 <- dy5 <-0
     }
   }

   ## this enables us specify constant food, instead of dynamic
   ##dy1<- v.f*dy1
   
   ##print(paste("animal is dead at time ",t,sep=""))
   
   list(c(dy1,dy2,dy3,dy4,dy5))
 }


PHI<-function(X.h, x){
  p<-x/(x+X.h)
  return(p)
}


food<-function(t, X.h, x){
  
  if(t>1) f<-(1-x/(x+X.h))*(2+sin(pi*t/2))##/2
  ##if(t%%7==0) f<-1000
  else f<-0
  
  f<-0
  return(f)
}

feed.switch<-function(t, y, p){
  c(-sin(2*pi*t/p["feed.space"]))
}

feed.map<-function(t, y, swID, p){
  if(swID==1){
    y[1]<-p["feed.level"]
  }
  return(y)
}


## This function takes data from the forward simulator (either full,
## or some subset), and applies the observation model and noise to
## generate "data" that would be like the observed data. The
## "observations" are then: L=alpha*y3 + eps.L (eps.L ~N(0, sds$L),
## and Negg=gamma*y5 + eps.L (eps.L ~N(0, sds$Negg), where y3 is
## length and y5 is reproduction.

## The name of this function should remain the same, since other
## functions elsewhere call it.

add.noise<-function(data, params, pred=FALSE){
  e.c<-0.0001
  ##alpha<-as.numeric(params["shape"])
  L.M<-as.numeric(params["L.M"])
  gamma<-as.numeric(params["gamma"])

  sd.L<-as.numeric(params["sd.L"])
  sd.Negg<-as.numeric(params["sd.Negg"])
  
  ##print(c(alpha, gamma))
  t<-data$time
  lt<-length(t)
  ##print(t)
  ##if(is.finite(alpha))  L<-data$y3*alpha
  ##else
  L<-data$y3*L.M

  Negg<-data$y5*gamma
  
  if(!pred){
    L<-rlnorm(lt, log(L+e.c), sd=sd.L) ##rnorm(t, data$y3*alpha, sd=sd.L)
    Negg<-rlnorm(lt, log(Negg+e.c), sd=sd.Negg) ##rnorm(t, data$y5*gamma, sd=sd.Negg)
  }
  ##w<-which(L<0)
  ##L[w]<-0
  ##w<-which(Negg<0)
  ##Negg[w]<-0
  
  return(list(t=t, L=L, Negg=Negg))

}


## Here are bits to set the parameters and initial conditions, and
## solve the system of odes, etc. There are also some plotting
## functions

## the parameter "var.food" indicates whether or not variable food
## levels are being considered -- 0 -> FALSE and 1-> TRUE. In other
## words, we can get a constant food environment by setting var.food=0

setparams.DEB<-function(J.EAm = 0.00361, y.EX = 0.5, v = 1.48,
                        k.M = 0.105, g = 0.498, kap = 0.1,
                        k.J.scale = 0.1, M.HP = 0.00575, 
                        gamma = 630, X.h = 0.16,
                        vol = 0.1, sd.L=0.1, sd.Negg=0.1,
                        l.M.HP=log(0.00575), F.eff= 0.2, 
                        L.M=4, var.food=1, feed.space=NULL, feed.level=NULL){
  
    params <- c(J.EAm = J.EAm, y.EX = y.EX, v = v, g = g,
                kap = kap, k.J.scale = k.J.scale, M.HP =
                M.HP,  L.M=L.M, gamma = gamma, X.h = X.h,
                vol = vol, sd.L=sd.L, sd.Negg=sd.Negg,
                l.M.HP=l.M.HP, F.eff=F.eff, 
                var.food=var.food, feed.space=feed.space,
                feed.level=feed.level)
    
  return(params)
}


## for this version of the model the initial enery density should be
## set to the "parental" density, i.e. to food(1, X.h, X0), but for
## the MCMC the "default is for e0==NULL, to indicate that e0 is
## effectively also "proposed by the algorithm".

setinits.DEB<-function(X0=0.2, e0=-Inf, L0 = 0.69, M.H0 = 0, M.R0 = 0){
  inits<-c(X0, e0, L0, M.H0, M.R0)
  return(inits)
}



## At this point this function has
## to be edited directly. Below are examples for how to build the list directly


##  all <- list()
##  all[[1]] <- list(params=c("a"), var=0.5, hyp=c(1,1), start=1)  ## 1-d proposal (if you give it a var command that is a vector, then a uniform proposal will be used instead of a normal)
##  all[[2]] <- list(params=c("b", "c"), var=diag(0.5, 2),
##                   hyp=matrix(c(1,2,3,4), nrow=2, byrow=TRUE),
##                   start=c(2,4)) ## joint 2-d proposal, does multivar norm proposals
##  all[[3]] <- list(params=c("d"), var=0) ## fix d


setup.mcmc.params<-function(hyper=NULL, w.p=NULL, samp.type=NULL,
                            prop.var=NULL, prop.mean=NULL, start=NULL,
                            joint=FALSE, n.joint=1, scale.var=NULL){


  ## first load the parameters for the priors
  if(is.null(hyper)) hyper<-make.hypers()
  
  ## Need to change these bits if parameters are only being proposed
  ## individually. Additionally, most of these can be changed by
  ## passing in a new vector with the appropriate name and format as
  ## an arguement to the function. However, most bits for any sets of
  ## parameters that are being proposed jointly must be changed
  ## manually below.
  
  if(is.null(w.p)) w.p <- c("sd.L", "sd.Negg", "kap", "g", "L.M", "v",
                            "l.M.HP") ##, "X.h")
  if(is.null(samp.type)) samp.type<-list(sd.L="ind", sd.Negg="ind", g="ind",
                                         k.M="ind", kap="ind", v="ind",
                                         l.M.HP="ind", shape="ind", X.h="ind",
                                         L.M="ind", F.eff="ind")
  if(is.null(prop.var)) prop.var<-list(sd.L=c(1,2), sd.Negg=c(1,2), g=0.75,
                                       k.M=0.45, kap=0.2, v=1.25,
                                       l.M.HP=5, shape=0.9, X.h=0.0002, L.M=0.1)
  if(is.null(prop.mean)) prop.mean<-list(sd.L=NULL, sd.Negg=NULL, g=NULL, k.M=NULL,
                                         kap=NULL, v=NULL, l.M.HP=NULL, shape=NULL,
                                         X.h=NULL, L.M=NULL)
  if(is.null(start)) start<-c(sd.L=0.128, sd.Negg=0.0718, l.M.HP=-5.285, shape=0.9,
                              k.M=0.822, kap=0.5, v=2.5, g=1.073, X.h=0.226, L.M=4,
                              F.eff=0.1)


  ## The following lines must be changed if parameters are being
  ## proposed jointly

  if(length(n.joint)==1) n.joint<-seq(1, n.joint, by=1)


  w1<-c("g", "v", "l.M.HP")
  cor1<-matrix(c( 1.0000000000000000000000,  0.8828502527611070682667,
-0.9462189953960029598079,  0.8828502527611070682667,
  1.0000000000000000000000, -0.8404684666239594648118,
 -0.9462189953960029598079, -0.8404684666239594648118,
  1.0000000000000000000000),
               ncol=3, byrow=TRUE)
  hyp1<-c(hyper$g, hyper$v, hyper$l.M.HP)
  m1<-c(0.5, 1.4, -5)
  t1<-"rw"

  w2<-c("L.M", "F.eff")
  cor2<-matrix(c(1,-0.9,
                 -0.9, 1), ncol=2, byrow=TRUE)
  hyp2<-c(hyper$L.M, hyper$F.eff)
  m2<-c(0.5, 2)
  t2<-"rw"

  w3<-c("kap", "v")
  cor3<-matrix(c(1,-0.8,
                 -0.8, 1), ncol=2, byrow=TRUE)
  hyp3<-c(hyper$g, hyper$v)
  m3<-c(0.5, 2)
  t3<-"rw"

  
  ws<-list(w1, w2, w3)
  ms<-list(m1, m2, m3)
  cors<-list(cor1, cor2, cor3)
  hyps<-list(hyp1, hyp2, hyp3)
  types<-list(t1, t2, t3)
  
  ## The code below here builds the appropriate structure with all the
  ## info needed by the mcmc to propose samples, etc.
  
  all<-list()
  
  if(joint){
    n.joint<-n.joint[!is.na(n.joint)]
    j<-0
    for(i in n.joint){
      j<-j+1
      all[[j]] <- list(params = ws[[i]],
                       ##var = diag(x=c(0.00055, 0.0005), nrow=2),
                       mean   = ms[[i]],
                       var    = make.sigNN(var=prop.var[ws[[i]]],
                         cor=cors[[i]],
                         scale=scale.var),
                       hyp    = matrix(hyps[[i]], nrow=length(ws[[i]]),
                         byrow=TRUE),
                       start  = start[ws[[i]]],
                       type   = types[[i]]
                       )
    }
  }
 
  if(!joint) n.joint<-NULL
  if(!is.null(w.p)){
    n.params <- length(n.joint)+length(w.p)
    
    for(i in (length(n.joint)+1):n.params){
      all[[i]]<-list(params=w.p[i-length(n.joint)])
    }
    
    for(i in 1:n.params){
      pp<-all[[i]]$params
      if(length(pp)==1){
        if(is.element(pp, w.p)){
          all[[i]]$type<-samp.type[[pp]]
          all[[i]]$mean<-prop.mean[[pp]]
          all[[i]]$var<-prop.var[[pp]]
          all[[i]]$hyp<-hyper[[pp]]
          all[[i]]$start<-start[[pp]]
        }
        else{
          all[[i]]$var<-0
          ##all[[i]]$hyp<-0
          ##all[[i]]$
        }
      }
    }
    
  }
  return(all)
}




## currently this only works to make a 2x2 cor matrix
make.sig22<-function(var=list(g=0.0009, eg=0.0009),
                   cor=-0.937){

  var<-as.numeric(unlist(var))

  cov<-cor*sqrt(var[1])*sqrt(var[2])
  
  sigma<-matrix(c(var[1], cov, cov, var[2]), ncol=2)

  return(sigma)
}

make.sigNN<-function(var=list(g=0.0009, eg=0.0009), cor, scale=NULL){

  if(!is.null(scale)) var.scale<-(2.38^2)/scale
  else var.scale<-1


  sd<-sqrt(as.numeric(unlist(var)))
  l<-length(var)
  
  sigma<-matrix(0, ncol=l, nrow=l)

  for(i in 1:l){
    for(j in 1:l){
      sigma[i,j]<-cor[i,j]*sd[i]*sd[j]*var.scale
    }
  }
  sigma<-as.matrix(nearPD(sigma)$mat)
  
  return(sigma)
}


