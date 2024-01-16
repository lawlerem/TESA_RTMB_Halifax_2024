library(RTMB)
df <- read.table("pol.dat", header=FALSE)
dat <- list(y=df[,1], X=as.matrix(df[,-1]))
par <- list(b=numeric(ncol(dat$X)), trans_a=0, logSigma=0, u=numeric(length(dat$y)))

f <- function(par){
  getAll(dat,par)
  timeSteps <- length(y)
  a <- 2*plogis(trans_a)-1
  sigma <- exp(logSigma)
  jnll <- -dnorm(u[1], 0, sigma/sqrt(1-a*a),log=TRUE)
  for(i in 2:timeSteps){
    jnll <- jnll - dnorm(u[i],a*u[i-1],sigma,log=TRUE)
  }
  eta <- X%*%b
  for(i in 1:timeSteps){
    lambda <- exp(eta[i]+u[i])   
    jnll <- jnll - dpois(y[i],lambda,log=TRUE) 
  }
  jnll
}                
            
obj <- MakeADFun(f,par,random="u")
opt<-nlminb(obj$par,obj$fn,obj$gr)

sdr<-sdreport(obj)
pl <- as.list(sdr,"Est")
plsd <- as.list(sdr,"Std")
