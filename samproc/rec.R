library(RTMB)
load("Robs.RData") # gets 'Robs' containing 'year', 'ssb', 'minAge' and 'Robs'
Robs$srmode <- 0 # later 1=ricker, 2=bh

par<-list(logSdR=0, logSdRobs=0, logR=numeric(length(Robs$year)))

f<-function(par){
  getAll(Robs, par)
  sdR <- exp(logSdR)
  sdO <- exp(logSdRobs)
  jnll <- 0
  for(i in 2:length(logR)){
    if(srmode==0){
      pred <- logR[i-1] 
    }else{
      stop(paste("srmode", srmode,"not implemented yet"))
    }
    jnll <- jnll - dnorm(logR[i], pred, sdR, log=TRUE)
  }
  for(i in 1:length(logR)){
    jnll <- jnll - dnorm(log(Robs[i]), logR[i], sdO, log=TRUE)
  }
  jnll
}

obj <- MakeADFun(f,par, random="logR", silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj)
pl <- as.list(sdr, "Est")
plsd <- as.list(sdr, "Std")

plot(Robs$year, log(Robs$Robs), xlab="Year", ylab="log(R)")
lines(Robs$year, pl$logR, col="darkred", lwd=3)
lines(Robs$year, pl$logR-2*plsd$logR, col="darkred", lwd=3, lty="dotted")
lines(Robs$year, pl$logR+2*plsd$logR, col="darkred", lwd=3, lty="dotted")
