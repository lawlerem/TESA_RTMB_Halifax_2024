library(RTMB)
dat <- list(Y=rnorm(1000,2,.3))
par <- list(mu=0, logsigma=0)
f <- function(par){
  getAll(par,dat)
  sigma <- exp(logsigma)
  nll <- -sum(dnorm(Y, mu, sigma, log=TRUE))
  ADREPORT(sigma)
  nll
}
obj <- MakeADFun(f, par)
opt <- nlminb(obj$par, obj$fn, obj$gr)
summary(sdreport(obj))
