library(RTMB)
load("allfleets.RData")

allfleets$keyQ <- rbind(c(NA,NA,NA,NA,NA,NA,NA,NA,NA),
                        c(NA, 1, 2, 3, 4, 5, 6, 7,NA),
                        c( 8, 9,10,11,12,13,NA,NA,NA))

allfleets$keySd <- rbind(c( 1, 1, 1, 1, 1, 1, 1, 1, 1),
                         c(NA, 2, 2, 2, 2, 2, 2, 2,NA),
                         c( 3, 3, 3, 3, 3, 3,NA,NA,NA))

allfleets$covType<-c(0,3,3)

allfleets$keyAR<-c(NA,1,2)

par <- list()
par$logQ <- numeric(max(allfleets$keyQ, na.rm=TRUE))
par$logsd <- numeric(max(allfleets$keySd, na.rm=TRUE))
par$tARphi <- numeric(sum(allfleets$covType==3))
par$missing <- numeric(sum(is.na(allfleets$obs)))


f<-function(par){
  getAll(allfleets,par)    
  nobs <- length(obs)
  sd <- exp(logsd)

  logobs <- log(obs)
  logobs[is.na(logobs)] <- missing

  nll <- 0
  
  logPred <- numeric(nobs)
  
  for(i in 1:nobs){
    y <- aux[i,1]-minYear+1
    f <- aux[i,2]
    a <- aux[i,3]-minAge+1
    Z <- F[y,a]+M[y,a]
    if(fleetTypes[f]==0){
      logPred[i] <- log(N[y,a])-log(Z)+log(1-exp(-Z))+log(F[y,a])
    }
    if(fleetTypes[f]==2){  
      logPred[i] <- logQ[keyQ[f,a]]+log(N[y,a])-Z*sampleTimes[f]
    }
    if(!(fleetTypes[f]%in%c(0,2))){  
      stop("This fleettype is has not been implemented yet")
    }
  }

  for(f in unique(aux[,2])){ # each fleet
    if(covType[f]==0){# independent      
      S <- diag(sd[na.omit(keySd[f,])]^2)
    }
    if(covType[f]==3){# AR
      S <- diag(sd[na.omit(keySd[f,])]^2)
      phi <- 2*plogis(tARphi[keyAR[f]])-1
      for(i in 2:ncol(S)){
        for(j in 1:(i-1)){
          S[i,j] <- sqrt(S[i,i])*sqrt(S[j,j])*(phi^(i-j))
          S[j,i] <- S[i,j]
        }
      }
    }
    if(!covType[f]%in%c(0,3)){#       
      stop("Covariance type not implemented")
    }

    
    #SU? res<-matrix(0,nrow=0,ncol=ncol(S)) 
    for(y in unique(aux[,1])){ # year within fleet 
      idx <- which((aux[,2]==f) & (aux[,1]==y))
      if(length(idx)!=0){
        nll <- nll - dmvnorm(logobs[idx],logPred[idx],S,log=TRUE)  
        #SU? res<-rbind(res,logobs[idx]-logPred[idx])
      }
    }
    #SU? nll<-nll-sum(dmvnorm(t(res),Sigma=S,log=TRUE))
  }
  REPORT(logPred)
  nll
}

obj <- MakeADFun(f, par, random="missing")
fit <- nlminb(obj$par, obj$fn, obj$gr)
est <- obj$report()$logPred

par(mfrow=c(1,3))
for(f in 1:3){
  idx<-which(allfleets$aux[,2]==f)
  matplot(xtabs(log(allfleets$obs[idx])~allfleets$aux[idx,1]
                +allfleets$aux[idx,3]), ylab="Log Obs")
  matplot(xtabs(est[idx]~allfleets$aux[idx,1]+allfleets$aux[idx,3]),
          type="l", add=TRUE)
}
