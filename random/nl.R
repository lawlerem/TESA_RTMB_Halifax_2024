library(RTMB)
obs <- as.matrix(read.table("nonlin.dat", header=FALSE))
par <- list(trans_phi=0, logsdLam=rep(0,2), logsdObs=rep(0,ncol(obs)), lam=matrix(0,nrow=nrow(obs),ncol=2))

f<-function(par){
  getAll(par)
  timeSteps <- nrow(lam)
  sdLam <- exp(logsdLam)
  sdObs <- exp(logsdObs)
  phi <- 2*plogis(trans_phi)-1

  jnll <- 0
  for(i in 2:timeSteps){    
    predU <- c(lam[i-1,1], phi*lam[i-1,2])
    jnll <- jnll - sum(dnorm(lam[i,],predU,sdLam,log=TRUE))
  }
  for(i in 1:timeSteps){
    predObs <- c(exp(lam[i,1])/(1+exp(lam[i,1])), lam[i,2], lam[i,1]*lam[i,2])
    jnll <- jnll -sum(dnorm(obs[i,],predObs,sdObs,log=TRUE))
  } 
  jnll    
}
obj <- MakeADFun(f,par,random="lam", silent=TRUE)
opt <- nlminb(obj$par,obj$fn,obj$gr)

sdr <- sdreport(obj)
pl <- as.list(sdr,"Est")
plsd <- as.list(sdr,"Std")
