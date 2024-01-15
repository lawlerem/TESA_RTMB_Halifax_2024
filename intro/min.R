library(RTMB)
par <- list(logk=c(0,0,0), logSigma=0)
dat <- read.table('min.dat', skip=3, head=FALSE)
nLogL <- function(par) {
  k <- exp(par$logk[1:3])
  sigma <- exp(par$logSigma)
  A <- matrix(c(-k[1], k[1], k[2], -(k[2] + k[3])),2,2); 
  x0 <- c(0, 100)
  sol <- function(t) 100 - sum(expm(A*t)%*%x0)
  pred <- numeric(nrow(dat))
  for(i in 1:nrow(dat))pred[i] <- sol(dat[i,1])
  -sum(dnorm(dat[, 2], mean = pred, sd = sigma, log = TRUE))
}

obj <- MakeADFun(nLogL,par, silent=T)

system.time(fit <- nlminb(obj$par,obj$fn, obj$gr))
