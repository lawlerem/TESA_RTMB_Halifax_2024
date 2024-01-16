library(RTMB)
dat <- list(y=scan("rwmissing.dat"))
par <- list(logSdRw=0, logSdObs=0, lam0=0, lam=numeric(length(dat$y)))

jnll<-function(par){
  getAll(par, dat)
  sdRw <- exp(logSdRw)
  sdObs <- exp(logSdObs)
  ret <- -dnorm(lam[1], mean=lam0, sd=sdRw, log=TRUE)
  ret <- ret - sum(dnorm(diff(lam), sd=sdRw, log=TRUE))
  idx <- which(!is.na(y))
  ret <- ret - sum(dnorm(y[idx], mean=lam[idx], sd=sdObs, log=TRUE))
  ret
}

obj <- MakeADFun(jnll, par, random="lam", silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj)
pl <- as.list(sdr, "Est")
plsd <- as.list(sdr, "Std")

plot(dat$y, xlab="Time", ylab="Y", las=1, ylim=range(pl$lam-2*plsd$lam,pl$lam+2*plsd$lam ))
lines(pl$lam, lwd=3, col="red")
lines(pl$lam-2*plsd$lam, lty="dotted", lwd=3, col="red")
lines(pl$lam+2*plsd$lam, lty="dotted", lwd=3, col="red")
