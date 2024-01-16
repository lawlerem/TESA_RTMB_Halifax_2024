library(RTMB)
dat <- list(y=scan("kf.dat"))
par <- list(logSdRw=0, logSdObs=0, lam0=0)

f<-function(par){
  getAll(par,dat)
    
  N <- length(y)  
  lam <- numeric(N)
  lamVar <- numeric(N)
  lamPred <- numeric(N) 
  lamPredVar <- numeric(N)
  lamSmooth <- numeric(N)
  lamSmoothVar <- numeric(N) 
  yPredVar <- numeric(N)
  w <- numeric(N)
  nll <- 0
  varLam <- exp(2.0*logSdRw)
  varY <- exp(2.0*logSdObs)
  lamPred[1] <- lam0
  lamPredVar[1] <- varLam
  yPredVar[1] <- lamPredVar[1]+varY
  w[1] <- y[1]-lamPred[1]
  lam[1] <- lamPred[1]+lamPredVar[1]/yPredVar[1]*w[1]
  lamVar[1] <- lamPredVar[1]-lamPredVar[1]/yPredVar[1]*lamPredVar[1]
  nll = nll-dnorm(w[1], sd=sqrt(yPredVar[1]),log=TRUE)

  for(i in 2:N){
    lamPred[i] <- lam[i-1]
    lamPredVar[i] <- lamVar[i-1]+varLam
    yPredVar[i] <- lamPredVar[i]+varY
    w[i] <- y[i]-lamPred[i]
    lam[i] <- lamPred[i]+lamPredVar[i]/yPredVar[i]*w[i]
    lamVar[i] <- lamPredVar[i]-lamPredVar[i]/yPredVar[i]*lamPredVar[i]
    nll <- nll-dnorm(w[i],sd=sqrt(yPredVar[i]), log=TRUE)
  }
  # ----------------- smoothing part 
  lamSmooth[N] <- lam[N]
  lamSmoothVar[N] <- lamVar[N]
  for(i in (N-1):1){
    varFrac <- lamVar[i]/lamPredVar[i+1]
    lamSmooth[i] <- lam[i]+varFrac*(lamSmooth[i+1]-lamPred[i+1])
    lamSmoothVar[i] <- lamVar[i]+varFrac*(lamSmoothVar[i+1]-lamPredVar[i+1])*varFrac
  }
  REPORT(lamSmooth)
  REPORT(lamSmoothVar)
  nll
}

obj<-MakeADFun(f,par)

obj$fn()
obj$gr()
opt<-nlminb(obj$par,obj$fn,obj$gr)

rep<-obj$report()

save(rep,file="kf.RData")
