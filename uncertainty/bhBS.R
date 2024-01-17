library(RTMB)
dat <- read.table("bh.dat", header=TRUE)
par <- list(loga=0, logb=0, logsigma=0)

nll<-function(par){
  getAll(par, dat)
  logR<-OBS(logR)  
  sigma <- exp(logsigma)
  pred <- loga+log(ssb)-log(1+exp(logb)*ssb)
  nll <- -sum(dnorm(logR,pred,sigma,TRUE))
  return(nll)
}

obj <- MakeADFun(nll, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
print(opt$par)

odat <- dat
doone <- function(){
  dat$logR <<- obj$simulate()$logR      
  objs <- MakeADFun(nll, par, silent=TRUE)
  opts <- nlminb(objs$par, objs$fn, objs$gr)
  opts$par
}
sim <- replicate(1000, doone())
dat <- odat

odat <- dat
doone<-function(){
  idx <- sample(1:nrow(dat), replace=TRUE)
  dat <<- odat[idx,]  
  objs <- MakeADFun(nll, par, silent=TRUE)
  opts <- nlminb(objs$par, objs$fn, objs$gr)
  opts$par
}
simbs<-replicate(1000, doone())
dat<-odat

par(mfrow=c(3,1))
for(i in 1:3){
  plot(density(sim[i,]), col="darkgreen", lwd=3)
  lines(density(simbs[i,]), col="orange", lwd=3)
}



