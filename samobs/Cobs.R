library(RTMB)
load("Cobs.RData")

par <- list(logsd=numeric(length(unique(Cobs$aux[,3]))))

f <- function(par){
  getAll(Cobs,par)
  nobs <- length(Cobs)
  sd <- exp(logsd)
  nll <- 0
  logPred <- numeric(nobs)
  
  for(i in 1:nobs){
    y <- aux[i,1]-minYear+1
    a <- aux[i,3]-minAge+1
    Z <- F[y,a]+M[y,a]
    logPred[i] <- log(N[y,a])-log(Z)+log(1-exp(-Z))+log(F[y,a])
    nll <- nll - dnorm(log(Cobs[i]),logPred[i],sd[a],log=TRUE)
  }
  REPORT(logPred)    
  nll
}

obj <- MakeADFun(f, par, map=list(logsd=factor(rep(1,length(par$logsd)))))
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr<-sdreport(obj)

matplot(rownames(Cobs$N), xtabs(log(Cobs$Cobs)~Cobs$aux[,1]+Cobs$aux[,3]), ylab="Log C", xlab="Year")
matplot(rownames(Cobs$N), xtabs(obj$report()$logPred~Cobs$aux[,1]+Cobs$aux[,3]), type="l", add=TRUE)

