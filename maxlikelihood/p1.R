library(RTMB)
dat <- list(X=2)
par <- list(p=.5)

f<-function(par){-dbinom(dat$X,100,par$p,log=TRUE)}

obj <- MakeADFun(f, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr, lower=c(0), upper=c(1))
summary(sdreport(obj))

#    Estimate Std. Error
# p      0.02 0.01398284

