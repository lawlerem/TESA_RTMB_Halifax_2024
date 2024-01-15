library(RTMB)
dat <- list()
dat$Y <- c(13, 5, 28, 28, 15, 4, 13, 4, 10, 17, 11, 13, 12, 17, 3)

par <- list()
par$logsize <- 0
par$logitp <- 0

nLogL <- function(par){
  -sum(dnbinom(dat$Y,exp(par$logsize),plogis(par$logitp),log=TRUE))
}

obj <- MakeADFun(nLogL, par)
opt <- nlminb(obj$par, obj$fn, obj$gr)
summary(sdreport(obj))
