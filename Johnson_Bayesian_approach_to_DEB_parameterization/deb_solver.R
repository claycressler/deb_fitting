## This file holds many of the functions needed to simulate the DEB
## model forward in time. It needs the package PBSddesolve to function.


require(PBSddesolve)

solve.DEB<-function(sim, params, inits, Tmax=400, numsteps=10000,
                    which="cont", sizestep=NULL, myswitch=NULL, mymap=NULL,
                    adj.inits=TRUE, adj.length=FALSE){
  
  if(adj.inits) inits<-adjust.inits(inits, params)

  if(is.null(sizestep) && is.null(numsteps)) stop("need to specify the times")
  if(is.null(sizestep)) times<- seq(0, Tmax, length=numsteps)
  else times<- seq(0, Tmax, by=sizestep)
  
  if(which=="cont"){## continuous feeding functions
    ##on.exit(freeglobaldata())
    out<-dde(y=inits, times=times, func=sim, parms=params)
  }
  else if(which=="discrete"){## discrete feeding situations
    ##on.exit(freeglobaldata())
    if(is.null(myswitch)|is.null(mymap)){
      print("not enough info to do the discrete feedings scenario")
      break
    }
    out<-dde(inits, times, sim, params, switchfunc=feed.switch, mapfunc=feed.map)
  }
  else stop("haven't chosen a correct solver scenario")

  if(adj.length){
    ## making sure to only return the output at the requested times
    out <- out[!duplicated(out$time),]
    out <- out[out$time %in% times,]
  }
  
  return(out)
}


adjust.inits<-function(inits, params){

  if(!is.finite(inits[1])) inits[1]<-as.numeric(params["F.eff"])
  
  if(!is.na(params["shape"])) inits[3]<-inits[3]/params["shape"]
  else inits[3]<-inits[3]/params["L.M"]
  
  if(!is.finite(inits[2])) inits[2]<-as.numeric(PHI(params["X.h"], inits[1]))
  return(inits)

}


test.sim <- function(sim, params, inits, ylim=c(0,600), Tmax=200,
                  numsteps=10000, which=2, scale=10){

  out<-solve.DEB(sim, params, inits, Tmax, numsteps, which)

  plot.DEB(out, scale)
  return(out)
 
}


plot.DEB<-function(out, scale=100, scaled.length=TRUE){
  par(mfrow=c(2,2))
   plot(out[,1],out[,2], type="l", lty=5, col="blue",
        lwd=2, xlab="time", ylab="Food")
   plot(out[,1],out[,3], type="l", col="red",  lty=4, lwd=2, xlab="time",
        ylab="scaled reserves")
  if(scaled.length){
    plot(out[,1],out[,4], type="l", col="green",lwd=2, xlab="time",
         ylab="scaled length")
  }
  else{
    plot(out[,1],out[,4], type="l", col="green",lwd=2, xlab="time",
         ylab="length")
  }
   plot(out[,1],out[,6]/scale, type="l", col=1, xlab="time",
        ylab=paste("maturity & reproduction/", scale, sep=""))
   lines(out[,1],out[,5], lty=2, col=2)
}


plot.DEB.red<-function(out){
  par(mfrow=c(2,1))
   plot(out$t,out$L, lty=5, col="blue",
        lwd=2, xlab="time", ylab="Length", pch=19)
   plot(out$t,out$Negg, col="red",  lty=4, lwd=2, xlab="time",
        ylab="Number of Eggs", pch=19)
}




## The following functions extract time series observations from the
## forward simulation. In particular, it extracts only those data
## points which are zero mod w.t

extract.data<-function(data, w.t=1, Tmax){

   if(length(w.t)==1){
    ww<-which(data$time%%w.t==0)
    ##print(ww)
    ##if(length(ww)!=(Tmax+1)){
    ##  print("something is wrong")
    ##  break
    ##}
    ##else
    data<-data[ww,]
  }
  else{
    ww<-NULL
    for(i in w.t){
      ww<-c(ww,which(data$time==i))
    }
    data<-data[ww,]
  }

  return(data)
  
}

## This function takes data from a forward simulation, extracts a
## subset of the data points using extract.data, and adds
## observational noise using add.noise

make.obs<-function(dt, params, w.t, Tmax){

  ##print(w.t)
  dt<-extract.data(dt, w.t, Tmax)
  #print(head(data))
  #return(data)
  dt<-add.noise(dt, params, pred=FALSE)
  return(dt)
  
}




## This function takes the simulator, params, inits, etc, and runs a
## forward simulation, and then extracts a subset of the points from
## the simulator (determined by w.t), without noise

make.states<-function(sim=DEB1, params, inits, Tmax, which=1, sizestep=0.01, w.t=1,
                      myswitch=NULL, mymap=NULL, numsteps=NULL){

  td<-td1<-0
  ##print(dt)

  ##if(!is.finite(inits[2])) inits[2]<-as.numeric(PHI(params["X.h"], inits[1]))
  
  td1<-solve.DEB(sim, params, inits, Tmax, numsteps=numsteps, which=which, sizestep,
                myswitch=myswitch, mymap=mymap)


  td<-extract.data(td1, w.t=w.t, Tmax)

  return(list(t=td$t, l=td$y3, n=td$y5))
  
}

make.pred<-function(sim=DEB1, params, inits, Tmax, which="const", sizestep=0.01, w.t=1,
                      myswitch=NULL, mymap=NULL, numsteps=NULL){

  temp<-temp1<-0
  ##print(dt)

  ##if(!is.finite(inits[2])) inits[2]<-as.numeric(PHI(params["X.h"], inits[1]))
  
  temp<-solve.DEB(sim, params, inits, Tmax, numsteps=numsteps, which=which, sizestep,
                myswitch=myswitch, mymap=mymap)

  ##print(w.t)
  temp1<-extract.data(temp, w.t, Tmax)
  #print(head(data))
  #return(data)
  temp<-add.noise(temp1, params, pred=TRUE)
  return(temp)
  
}


##################################################################
## These functions are no longer being used
##################################################################

discrete.feed<-function(inits, times, sim, params, feed.space, feed.level){
  out<-temp<-0
  t.m<-max(times)
  ww<-which(times%%feed.space==0)
  l<-length(ww)
  ##print(c(t.m, times[ww[l]]))
  temp<-dde(y=inits, times=times[1:ww[2]], func=sim, parms=params)

  d<-dim(temp)
  out<-data.frame(matrix(0,nrow=length(times), ncol=d[2]))

  names(out)<-names(temp)
  out[,1]<-times
  out[1:ww[2],2:d[2]]<-temp[,2:d[2]]

  freeglobaldata()
  
  for(i in 3:l){

    temp.inits<-as.numeric(out[ww[i-1],2:d[2]])

    temp.inits[1]<-feed.level
    
    temp<-dde(y=temp.inits, times=times[(ww[i-1]+1):ww[i]], func=sim, parms=params)
    
    ##print(c(times[(ww[i-1]+1)], times[ww[i]]))
    ##print(head(temp))
    ##print(tail(temp))

    out[(ww[i-1]+1):ww[i],2:d[2]]<-temp[,2:d[2]]
    freeglobaldata()

  }
  if(t.m>times[ww[l]]){

    temp.inits<-as.numeric(out[ww[l],2:d[2]])
    temp.inits[1]<-feed.level
    temp<-dde(y=temp.inits, times=times[(ww[l]+1):length(times)], func=sim, parms=params)
    
    ##print(c(times[(ww[l]+1)], times[length(times)]))
    ##print(head(temp))
    ##print(tail(temp))

    out[(ww[l]+1):length(times),2:d[2]]<-temp[,2:d[2]]
    freeglobaldata()

  }

  temp<-0
  
  return(out)

}

