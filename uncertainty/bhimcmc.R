library(RTMB)
dat <- read.table("bh.dat", header=TRUE)
dat0 <- read.table("bh0.dat", header=TRUE)
dat<-rbind(dat,dat0)
par <- list(loga=0, logb=0, logsigma=0)

nll<-function(par){
  getAll(par, dat)
  sigma <- exp(logsigma)
  pred <- loga+log(ssb)-log(1+exp(logb)*ssb)
  nll <- -sum(dnorm(logR,pred,sigma,TRUE))
  return(nll)
}

obj <- RTMB::MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- summary(sdreport(obj))

library(tmbstan)
fitmcmc2 <- tmbstan(obj, chains=1,
              iter=10000, init=list(opt$par),
              lower=c(-50,-50, -50), upper=c(50,50,50) )

mc <- extract(fitmcmc2, pars=names(obj$par),
              inc_warmup=TRUE, permuted=FALSE)

#npar<-dim(mc)[3]
#layout(rbind(rep(1,npar),c(2:(npar+1)))) 
#matplot(mc[,1,], type="l")
#
#for(i in 1:npar){
#  pn <- attr(mc,"dimnames")$parameters[i]
#  plot(density(mc[,1,i]), main=pn, col=i, lwd=3)
#  abline(v=sdr[i,1], lwd=3)
#  xx <- pretty(range(mc[,1,i]),100)
#  lines(xx,dnorm(xx,sdr[i,1],sdr[i,2]),
#        col=gray(.5), lwd=3)
#}

mci<-mc[,1,]
