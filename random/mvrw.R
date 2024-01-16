library(RTMB)
obs <- read.table("files/mvrw.dat", header=FALSE)
par <- list(logSdLam=rep(0,ncol(obs)), logSdObs=rep(0,ncol(obs)), lambda=matrix(0, nrow=nrow(obs), ncol=ncol(obs)))

f<-function(par){
  getAll(par)
  stateDim <- ncol(lambda)
  timeSteps <- nrow(lambda)
  sdLam <- exp(logSdLam)
  sdObs <- exp(logSdObs);
  S <- diag(sdLam^2) 
  Q <- diag(sdObs^2) 
  ret <- 0
  for(i in 2:timeSteps){
    ret <- ret - dmvnorm(lambda[i,], lambda[i-1,], Sigma=S, log=TRUE)
  }
  for(i in 1:timeSteps){
    ret <- ret - dmvnorm(obs[i,], lambda[i,], Sigma=Q, log=TRUE)
  } 
  ret
}

obj <- MakeADFun(f,par,random="lambda")
opt<-nlminb(obj$par,obj$fn,obj$gr)

