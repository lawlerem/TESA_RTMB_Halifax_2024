library(RTMB)
obs <- read.table("files/mvrw.dat", header=FALSE)
par <- list(logSdLam=rep(0,ncol(obs)), logSdObs=rep(0,ncol(obs)), lambda=matrix(0, nrow=nrow(obs), ncol=ncol(obs)), tRho=0)

f <- function(par){
  getAll(par)
  stateDim <- ncol(lambda)
  timeSteps <- nrow(lambda)
  sdLam <- exp(logSdLam)
  sdObs <- exp(logSdObs)
  rho <- 2*plogis(tRho)-1
  S <- diag(sdLam^2)
  for(i in 1:(stateDim-1)){
    for(j in (i+1):stateDim){
      S[i,j] <- rho^(j-i)*sdLam[i]*sdLam[j]
      S[j,i] <- S[i,j]  
    }
  }
  Q <- diag(sdObs^2) 
  ret <- -sum(dmvnorm(lambda[-1,], lambda[-nrow(lambda),], Sigma=S, log=TRUE))
  ret <- ret - sum(dmvnorm(obs, lambda, Sigma=Q, log=TRUE))
  ret
}
obj <- MakeADFun(f,par,random="lambda", silent=TRUE)
opt <- nlminb(obj$par,obj$fn,obj$gr)
sdr <- sdreport(obj)
pl <- as.list(sdr,"Est")
plsd <- as.list(sdr,"Std")

Time <- 1:nrow(obs)
matplot(Time, obs, type="n")
for(j in 1:ncol(obs))polygon(cbind(c(Time,rev(Time)), c(pl$lambda[,j]-2*plsd$lambda[,j], rev(pl$lambda[,j]+2*plsd$lambda[,j]))),
                             col=scales:::alpha(palette()[j],.2), border=NA)
matplot(pl$lambda, type="l", lwd=2, add=TRUE, lty="solid")
matplot(Time, obs, pch=4, add=TRUE)
