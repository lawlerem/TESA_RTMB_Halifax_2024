library(RTMB)
dat <- list()
dat$Y <- c(13, 5, 28, 28, 15, 4, 13, 4, 10, 17, 11, 13, 12, 17, 3) 

par <- list()
par$logsize <- 0
par$logitp <- 0

nll<-function(par){
  getAll(par)
  size <- exp(logsize)
  p <- plogis(logitp)
  -sum(dnbinom(dat$Y, size=size, prob=p, log=TRUE))
}

obj <- RTMB::MakeADFun(nll, par)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- summary(sdreport(obj))

library(tmbstan)
fitmcmc2 <- tmbstan(obj, chains=1,
              iter=10000, init=list(opt$par),
              lower=c(-5,-5), upper=c(5,5) )

mc <- extract(fitmcmc2, pars=names(obj$par),
              inc_warmup=TRUE, permuted=FALSE)

npar<-dim(mc)[3]
layout(rbind(rep(1,npar),c(2:(npar+1)))) 
matplot(mc[,1,], type="l")

for(i in 1:npar){
  pn <- attr(mc,"dimnames")$parameters[i]
  plot(density(mc[,1,i]), main=pn, col=i, lwd=3)
  abline(v=sdr[i,1], lwd=3)
  xx <- pretty(range(mc[,1,i]),100)
  lines(xx,dnorm(xx,sdr[i,1],sdr[i,2]),
        col=gray(.5), lwd=3)
}


