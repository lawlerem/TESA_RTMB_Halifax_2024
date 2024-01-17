library(RTMB)
x<-scan("ar1.dat", quiet=TRUE)

dat <- list(x=x, code=0)
par <- list(logSigma=0, tPhi=1)

f<-function(par){
  getAll(dat,par)
  phi <- 2*plogis(tPhi)-1
  timeSteps <-length(x)
  sd <- exp(logSigma)

  nll=0;

  if(code==0){
    nll <- nll -dnorm(x[1],0,sqrt(sd*sd/(1-phi*phi)),log=TRUE)
    for(i in 2:timeSteps){    
      nll <- nll -dnorm(x[i],phi*x[i-1],sd,log=TRUE)
    }
  }
  
  if(code==1){
    nll <- nll - dautoreg(x,phi=phi,scale=sqrt(sd*sd/(1-phi*phi)), log=TRUE)
  }

  if(code==2){
    S <- matrix(0,timeSteps,timeSteps)
    for(i in 1:timeSteps){
      for(j in 1:timeSteps){
        S[i,j] <- sd*sd/(1-phi*phi)*phi^abs(i-j)
      }
    }
    nll <- nll - dmvnorm(x,0,S,log=TRUE)
  }

  if(code==3){
    Q<-Matrix::spMatrix(timeSteps,timeSteps)
    diag(Q)<-c(1,rep(1+phi^2, timeSteps-2),1)
    Q[cbind(2:timeSteps,1:(timeSteps-1))] <- -phi
    Q[cbind(1:(timeSteps-1), 2:timeSteps)] <- -phi
    Q <- Q/sd/sd
    nll <- nll - dgmrf(x, Q=Q, log=TRUE)
  }
  ADREPORT(phi) 
  nll
}

dat$code=0
obj <- MakeADFun(f,par, silent=TRUE)
system.time(opt<-nlminb(obj$par,obj$fn,obj$gr))
c(opt$obj, opt$par)
dat$code=1
obj <- MakeADFun(f,par, silent=TRUE)
system.time(opt<-nlminb(obj$par,obj$fn,obj$gr))
c(opt$obj, opt$par)
dat$code=2
obj <- MakeADFun(f,par, silent=TRUE)
system.time(opt<-nlminb(obj$par,obj$fn,obj$gr))
c(opt$obj, opt$par)
dat$code=3
obj <- MakeADFun(f,par, silent=TRUE)
system.time(opt<-nlminb(obj$par,obj$fn,obj$gr))
c(opt$obj, opt$par)

