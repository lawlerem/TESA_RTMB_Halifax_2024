library(RTMB)
load("Nobs.RData")
Nobs$srmode <- 0
matplot(log(Nobs$Nobs), main="logN")

par <- list()
par$logsdR <- 0
par$logsdS <- 0
par$logsdO <- 0
par$logN <- matrix(0, nrow=nrow(Nobs$Nobs), ncol=ncol(Nobs$Nobs))

f<-function(par){
  getAll(par, Nobs)
  nrow <- nrow(Nobs)
  ncol <- ncol(Nobs)
  sdR <- exp(logsdR)
  sdS <- exp(logsdS)
  sdO <- exp(logsdO)  

  jnll <- 0 
  
  for(y in 2:nrow){
    if(srmode==0){
      pred <- logN[y-1,1]
    }else{
      stop(paste("srmode", srmode,"not implemented yet"))
    }      
    jnll <- jnll - dnorm(logN[y,1],pred,sdR,log=TRUE)
  }

  for(y in 2:nrow){
    for(a in 2:ncol){
      pred <- logN[y-1,a-1]-F[y-1,a-1]-M[y-1,a-1]
      if(a==ncol){
        pred <- log(exp(pred)+exp(logN[y-1,a]-F[y-1,a]-M[y-1,a]))
      }
      jnll <- jnll - dnorm(logN[y,a],pred,sdS,log=TRUE)
    }
  }

  for(y in 1:nrow){
    for(a in 1:ncol){
      jnll <- jnll - dnorm(log(Nobs[y,a]),logN[y,a],sdO,log=TRUE)
    }
  }
  
  jnll
}

obj <- MakeADFun(f, par, random="logN", silent=TRUE)

fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
pl <- as.list(sdr,"Est")
matplot(pl$logN, type="l", add=TRUE)


