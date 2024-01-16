load("fsa.RData") # gets "dat"

set.seed(1234)
dat$obs[sample(1:400, 10)]<-100000000

drob <- function(x,m,s, p=0.05, log=FALSE){
  z<-(x-m)/s  
  res<-((1-p)*dnorm(z)+p*dt(z,df=3))/s
  if(log){
    return(log(res))
  }else{
    return(res)
  }    
}

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
  
  ans <- -sum(drob(logObs,logPred,sdvec,log=TRUE))

  logssb <- log(apply(exp(logN)*stockMeanWeight*propMature,2,sum))
  ADREPORT(logssb)
  REPORT(F)
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




