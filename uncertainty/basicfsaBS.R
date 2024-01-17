load("fsa.RData") # gets "dat"

library(RTMB)

par <- list(
  logN1Y=rep(0,nrow(dat$M)),
  logN1A=rep(0,ncol(dat$M)-1),
  logFY=rep(0,ncol(dat$M)),
  logFA=rep(0,nrow(dat$M)),
  logSdCatch=0, 
  logQ=rep(0,length(unique(dat$age[dat$fleet==2]))),
  logSdSurvey=0
)

nll<-function(par){
  getAll(par, dat)
  obs<-OBS(obs)
  na <- max(age)-min(age)+1
  ny <- max(year)-min(year)+1

  ## setup F
  F <- outer(exp(logFA),exp(logFY))

  ## setup N
  logN <- matrix(0, nrow=na, ncol=ny)
  logN[,1] <- logN1Y  
  for(y in 2:ny){
    logN[1,y] <- logN1A[y-1]
    for(a in 2:na){
      logN[a,y] <- logN[a-1,y-1]-F[a-1,y-1]-M[a-1,y-1]
    }
  } 

  # Match to observations
  logObs <- log(obs)
  logPred <- numeric(length(logObs))
  sdvec <- numeric(length(logObs))
  for(i in 1:length(logObs)){
    a <- age[i]-min(age)+1
    y <- year[i]-min(year)+1
    if(fleet[i]==1){
      logPred[i] <- log(F[a,y])-log(F[a,y]+M[a,y])+log(1-exp(-F[a,y]-M[a,y]))+logN[a,y]
      sdvec[i] <- exp(logSdCatch)
    }else{
      logPred[i] <- logQ[a]-(F[a,y]+M[a,y])*surveyTime+logN[a,y]
      sdvec[i] <- exp(logSdSurvey)
    }    
  }
  
  ans <- -sum(dnorm(logObs,logPred,sdvec,TRUE))
  
  logssb <- log(apply(exp(logN)*stockMeanWeight*propMature,2,sum))
  ADREPORT(logssb)
  REPORT(logPred)
  REPORT(logObs)  
  return(ans)
}

obj <- MakeADFun(nll, par, map=list(logFA=factor(c(1:4,NA,NA,NA))), silent=TRUE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(iter.max=1000,eval.max=1000))
sdrep <- sdreport(obj)
pl <- as.list(sdrep, "Est", report=TRUE)
plsd <- as.list(sdrep, "Std", report=TRUE)

yr<-sort(unique(dat$year))
plot(yr, exp(pl$logssb), type="l", lwd=5, col="red", ylim=c(0,550000), xlab="Year", ylab="SSB")
lines(yr, exp(pl$logssb-2*plsd$logssb), type="l", lwd=1, col="red")
lines(yr, exp(pl$logssb+2*plsd$logssb), type="l", lwd=1, col="red")

##  odat<-dat
##  ssbsimBS<-list()
##  for(i in 1:100){
##    idx<-sample(1:length(odat$obs), replace=TRUE)
##    dat$obs<-odat$obs[idx]
##    dat$fleet<-odat$fleet[idx]
##    dat$year<-odat$year[idx]
##    dat$age<-odat$age[idx]  
##    objnew <- MakeADFun(nll, par, map=list(logFA=factor(c(1:4,NA,NA,NA))), silent=TRUE)
##    opt <- nlminb(objnew$par, objnew$fn, objnew$gr, control=list(iter.max=1000,eval.max=1000))
##    sdrepnew <- sdreport(objnew)
##    plnew <- as.list(sdrepnew, "Est", report=TRUE)
##    ssbsimBS[[i]]<-exp(plnew$logssb)
##    lines(yr, exp(plnew$logssb), col=gray(.2), lwd=.5)
##  }
##  dat<-odat

odat<-dat
i1<-which(odat$fleet==1)
i2<-which(odat$fleet==2)
p<-obj$report()$logPred
o<-obj$report()$logObs
r<-o-p
ssbsimBS<-list()
for(i in 1:100){
  dat$obs[i1]<-exp(p[i1]+sample(r[i1], replace=TRUE))
  dat$obs[i2]<-exp(p[i2]+sample(r[i2], replace=TRUE))
  objnew <- MakeADFun(nll, par, map=list(logFA=factor(c(1:4,NA,NA,NA))), silent=TRUE)
  opt <- nlminb(objnew$par, objnew$fn, objnew$gr, control=list(iter.max=1000,eval.max=1000))
  sdrepnew <- sdreport(objnew)
  plnew <- as.list(sdrepnew, "Est", report=TRUE)
  ssbsimBS[[i]]<-exp(plnew$logssb)
  lines(yr, exp(plnew$logssb), col=gray(.2), lwd=.5)
}
dat<-odat

lines(yr, exp(pl$logssb-2*plsd$logssb), type="l", lwd=1, col="red")
lines(yr, exp(pl$logssb+2*plsd$logssb), type="l", lwd=1, col="red")
matplot(yr, t(apply(do.call(rbind, ssbsimBS), 2, quantile, prob=c(0.025, 0.975))), add=TRUE, col="orange", type="l", lwd=2, lty="solid")
 
