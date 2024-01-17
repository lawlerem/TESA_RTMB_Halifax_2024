library(RTMB)
load("Nobs2.RData")
Nobs$srmode <- 2
matplot(log(Nobs$Nobs), main="logN")

par <- list()
par$logsdR <- 0
par$logsdS <- 0
par$logsdO <- 0
par$logN <- matrix(0, nrow=nrow(Nobs$Nobs), ncol=ncol(Nobs$Nobs))
par$rickerpar <- if(Nobs$srmode==1){c(1,-10)}else{numeric(0)}
par$bhpar <- if(Nobs$srmode==2){c(1,1)}else{numeric(0)}


ssbFUN <- function(logN, F, M, SW, MO, PF, PM){
  nrow <- nrow(logN)
  ncol <- ncol(logN)
  ret <- numeric(nrow)
  for(y in 1:nrow){
    for(a in 1:ncol){
      ret[y] = ret[y]+SW[y,a]*MO[y,a]*exp(logN[y,a])*exp(-PF[y,a]*F[y,a]-PM[y,a]*M[y,a])
    }
  }
  return(ret);
}


f<-function(par){
  getAll(par, Nobs)
  nrow <- nrow(Nobs)
  ncol <- ncol(Nobs)
  sdR <- exp(logsdR)
  sdS <- exp(logsdS)
  sdO <- exp(logsdO)  
  minAge <- min(age)
  ssb <- ssbFUN(logN,F,M,SW,MO,PF,PM)
  
  jnll <- 0 
  
  for(y in 2:nrow){
    thisSSB <- ifelse((y-minAge-1)>(-.5),ssb[y-minAge],ssb[1])
      
    if(srmode==0){ #RW
      pred <- logN[y-1,1]
    }
    if(srmode==1){ #Ricker
      pred <- rickerpar[1]+log(thisSSB)-exp(rickerpar[2])*thisSSB
    }
    if(srmode==2){ #BH
      pred <- bhpar[1]+log(thisSSB)-log(1.0+exp(bhpar[2])*thisSSB)
    }

    if(!(srmode%in%c(0,1,2))){
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

obj <- MakeADFun(f, par, random="logN", silent=FALSE)

fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
pl <- as.list(sdr,"Est")
matplot(pl$logN, type="l", add=TRUE)


