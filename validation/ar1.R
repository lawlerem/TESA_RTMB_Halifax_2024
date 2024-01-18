library(RTMB)
load("cpue.RData")
dat = list(y = y)
par = list(mu=0, logSigma = 0, phiTrans = 0, gamma = rep(0,length(dat$y)))

f<-function(par){
  getAll(dat,par)
  phi <- 2*plogis(phiTrans)-1
  sd <- exp(logSigma) 
  jnll <- -dautoreg(gamma, mu=mu, phi=phi, scale=sd, log=TRUE)   
  jnll <- jnll -sum(dpois(y,exp(gamma),log=TRUE))
  jnll
}

obj <- MakeADFun(f,par,random = "gamma")

fit <- nlminb(obj$par,obj$fn,obj$gr)
