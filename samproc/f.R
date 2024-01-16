library(RTMB)
load("Fobs.RData")
matplot(Fobs$year, log(Fobs$Fobs), xlab="Year", ylab="logF", pch=colnames(Fobs$Fobs))

Fobs$fcormode <- 2

par <- list()
par$logsdF <- rep(0,ncol(Fobs$Fobs))
par$transRhoF <- if(Fobs$fcormode==0){numeric(0)}else{.1}
par$logsdO <- 0
par$logF <- matrix(0, nrow=nrow(Fobs$Fobs), ncol=ncol(Fobs$Fobs))

f <- function(par){
  getAll(Fobs,par)
    
  nrow <- nrow(Fobs)
  ncol <- ncol(Fobs)
  
  sdF <- exp(logsdF)
  sdO <- exp(logsdO)  
  rho <- 0
  
  jnll <- 0
  
  SigmaF <- matrix(0, ncol, ncol)

  if(fcormode==0){
    diag(SigmaF) <- sdF*sdF
  }

  if(fcormode==1){
    diag(SigmaF) <- sdF*sdF
    rhoF <- 2*plogis(transRhoF[1])-1
    for(i in 2:ncol){
      for(j in 1:(i-1)){
        SigmaF[i,j] <- rhoF*sdF[i]*sdF[j]
        SigmaF[j,i] <- SigmaF[i,j] 
      }
    }
  }
  
  if(fcormode==2){
    diag(SigmaF) <- sdF*sdF
    rhoF <- 2*plogis(transRhoF[1])-1
    for(i in 2:ncol){
      for(j in 1:(i-1)){
        SigmaF[i,j] <- sdF[i]*sdF[j]*(rhoF^(i-j))
        SigmaF[j,i] <- SigmaF[i,j] 
      }
    }
  }

  for(y in 2:nrow){
    jnll <- jnll - dmvnorm(logF[y,],logF[y-1,],SigmaF,log=TRUE)
  }

  for(y in 1:nrow){
    for(a in 1:ncol){
      jnll <- jnll - dnorm(log(Fobs[y,a]), logF[y,a], sdO, log=TRUE)
    }
  }
  jnll    
}

obj <- MakeADFun(f, par, random="logF", map=list(logsdF=factor(rep(1,ncol(Fobs$Fobs)))))
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdr<-sdreport(obj)
matplot(Fobs$year, as.list(sdr,"Est")$logF, type="l", add=TRUE)

## res
## ID: nll=441.9986
## CS: nll=278.1695
## AR: nll=230.7495
## parallel: nll ~ 485




#/* For the AR(1) model it is possible to write the process part as: 
# 
#    using namespace density;
#    VECSCALE_t<AR1_t<N01<Type> > > nlAR=VECSCALE(AR1(itrans(transRho(0))),sdF);
#    for(int y=1; y<nrow; ++y){
#      jnll += nlAR(logF.row(y)-logF.row(y-1));
#    }
#
#*/
