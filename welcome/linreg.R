library(RTMB)
dat <- read.table("linreg.dat", header=TRUE)

nll<-function(par){
  getAll(dat,par)
  pred <- alpha+beta*x
  -sum(dnorm(y,pred,exp(logSigma),TRUE))
}

par <- list(alpha=0, beta=0, logSigma=0)

obj <- MakeADFun(nll,par)
opt <- nlminb(obj$par,obj$fn,obj$gr)
sdrep <- sdreport(obj)
summary(sdrep)

