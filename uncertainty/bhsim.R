library(RTMB)
dat <- read.table("bh.dat", header=TRUE)
par <- list(loga=0, logb=0, logsigma=0)

nll<-function(par){
  getAll(par, dat)
  logR<-OBS(logR)  
  sigma <- exp(logsigma)
  pred <- loga+log(ssb)-log(1+exp(logb)*ssb)
  REPORT(pred)
  nll <- -sum(dnorm(logR,pred,sigma,TRUE))
  return(nll)
}

obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)

plot(dat$ssb, obj$report()$pred, type="l")

print(opt$par)

oopt<-opt
odat<-dat
sim<-list()
simp<-list()
for(i in 1:100){
  dat$logR<-obj$simulate()$logR      
  newobj <- MakeADFun(nll, par, silent=TRUE)
  newopt<-nlminb(newobj$par, newobj$fn, newobj$gr)
  sim[[i]]<-newopt$par
  simp[[i]]<-newobj$report()$pred
}

invisible(lapply(simp, function(y)lines(dat$ssb, y, col="gray")))

boxplot(do.call(rbind, sim))
points(1:3, oopt$par, pch=4, lwd=3, cex=3, col="hotpink")







