library(RTMB)
dat <- list(Y=c(13,5,28,28,15,4,13,4,10,17,11,13,12,17,3)) 
par <- list(logsize=0, logitp=0)

nll<-function(par){
  getAll(par)
  size <- exp(logsize)
  p <- plogis(logitp)
  -sum(dnbinom(dat$Y, size=size, prob=p, log=TRUE))
}

obj <- RTMB::MakeADFun(nll, par)
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)
sdr <- summary(rep)

est <- opt$par
npar <- length(est)
cov <- rep$cov.fixed # or just # diag(npar) #  

nosim <- 10000
mcmc <- matrix(NA,nrow=nosim, ncol=npar)
mcmc[1,] <- est
colnames(mcmc) <- names(est)

steps <- MASS::mvrnorm(nosim, mu=rep(0,npar), Sigma=cov)
U <- runif(nosim)

currentValue <- obj$fn(mcmc[1,])
for(i in 2:nosim){
  proposal <- mcmc[i-1,]+steps[i,]
  proposalValue <- obj$fn(proposal)
  if(U[i]<exp(currentValue-proposalValue)){
    mcmc[i,] <- proposal
    currentValue <- proposalValue
  }else{
    mcmc[i,] <- mcmc[i-1,]
  }
}

layout(rbind(rep(1,npar),c(2:(npar+1)))) 
matplot(mcmc, type="l")
for(i in 1:npar){
  pn <- names(est)[i]
  plot(density(mcmc[,i]), main=pn, col=i, lwd=3)
  abline(v=sdr[i,1], lwd=3)
  xx <- pretty(range(mcmc[,i]),100)
  lines(xx,dnorm(xx, sdr[i,1],sdr[i,2]), col=gray(.5), lwd=3)
}
