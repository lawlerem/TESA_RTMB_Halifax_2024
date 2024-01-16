library(RTMB)
z1 <- rnorm(1000)
z2 <- rnorm(1000)
Z <- rbind(z1,z2)
A <- rbind(c(1,0),c(-1,1))
b <- c(4,-2)
X <- A%*%Z+b

dat <- list(X=t(X))
par <- list()
par$mu <- c(0,0)
par$logSigma <- c(0,0)
par$tRho <- 1

f <- function(par){
  getAll(par,dat)  
  sigma <- exp(logSigma)
  rho <- 2*plogis(tRho)-1
  Sigma <- rbind( c( sigma[1]^2,             sigma[1]*sigma[2]*rho),
                  c( sigma[1]*sigma[2]*rho,  sigma[2]^2)            )
  REPORT(Sigma)
  -sum(dmvnorm(X, mu, Sigma=Sigma, log=TRUE))
}

obj <- MakeADFun(f, par)
opt <- nlminb(obj$par, obj$fn, obj$gr)
summary(sdreport(obj))


