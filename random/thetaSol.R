library(RTMB)
Y<-scan('theta.dat', quiet=TRUE)
par <- list(X=numeric(length(Y)), logr0=0, logtheta=0, logK=6, logQ=0, logR=0)

f<-function(par){
  getAll(par)
  r0 <- exp(logr0)
  theta <- exp(logtheta)
  K <- exp(logK)
  Q <- exp(logQ)
  R <- exp(logR)
  timeSteps <- length(Y)
  jnll <- 0
  for(i in 2:timeSteps){
    pred <- X[i-1]+r0*(1.0-(exp(X[i-1])/K)^theta)
    jnll <- jnll - dnorm(X[i],pred,sqrt(Q),log=TRUE)
  }
  jnll <- jnll - sum(dnorm(Y,X,sqrt(R),log=TRUE))
  jnll
}
obj <- MakeADFun(f,par,random="X", silent=TRUE)
opt <- nlminb(obj$par,obj$fn,obj$gr)

sdr<-sdreport(obj)
pl <- as.list(sdr,"Est")
plsd <- as.list(sdr,"Std")
plot(Y); lines(pl$X, col="hotpink", lwd=3)
lines(pl$X-2*plsd$X, col="hotpink", lwd=3, lty="dotted"); lines(pl$X+2*plsd$X, col="hotpink", lwd=3, lty="dotted")
