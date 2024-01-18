library(RTMB)
load("cpue.RData")
dat <- list(y = y) # rm(y)
par <- list(mu=0, logSigma = 0, phiTrans = 0, gamma = rep(0,length(dat$y)))

f<-function(par){
  getAll(dat,par)
  y<-OBS(y)  
  phi <- 2*plogis(phiTrans)-1
  sd <- exp(logSigma) 
  jnll <- -dautoreg(gamma, mu=mu, phi=phi, scale=sd, log=TRUE)   
  jnll <- jnll -sum(dpois(y,exp(gamma),log=TRUE))
  REPORT(mu)
  REPORT(phi)
  stepsd <- sqrt(sd^2*(1-phi^2)) # because sd=sqrt(stepsd^2*(1-phi^2))  
  REPORT(stepsd)  
  jnll
}

obj <- MakeADFun(f,par,random = "gamma")

fit <- nlminb(obj$par,obj$fn,obj$gr)

par(mfrow=c(1,2))

# RES1
res <- oneStepPredict(obj,method="oneStepGeneric", discrete=TRUE, range=c(0,Inf))
qqnorm(res$residual)
abline(0,1)

# RES2

sdr <- sdreport(obj)
estX <- summary(sdr,"random")
C <- solve(obj$env$spHess(obj$env$last.par.best , random=TRUE))
gamma.star <- MASS::mvrnorm (1, estX[,1] ,C)

rep<-obj$report()

#(gamma[i]-mu) = phi*(gamma[i-1]-mu) + normal dist with sd=stepsd

z<-((gamma.star[-1]-rep$mu)-rep$phi*(gamma.star[-length(gamma.star)]-rep$mu))/rep$stepsd
qqnorm(z)
abline(0,1)

# Jit

doone<-function()nlminb(obj$par+rnorm(length(obj$par), sd=.1),obj$fn,obj$gr)$par
jit<- replicate(100, doone())
boxplot(t(jit))

# sim study full
odat <- dat
doone<-function(){
    dat <<- list()
    dat$y <<- obj$simulate()$y
    objsim <- MakeADFun(f,par,random = "gamma", silent=TRUE)
    fitsim <- nlminb(objsim$par,objsim$fn,objsim$gr)
    fitsim$par
}
sim <- replicate(100, doone())
dat <- odat
boxplot(t(sim))
points(1:length(fit$par), fit$par, cex=5, pch=4, lwd=3, col="darkgreen")

# sim study full conditional

# sim study full
odat <- dat
doone<-function(){
    dat <- list()
    sdr <- sdreport(obj)
    pl <- as.list(sdr, "Est") 
    objcon <- MakeADFun(f,pl, silent=TRUE)
    dat$y <<- objcon$simulate()$y
    objsim <- MakeADFun(f,par,random = "gamma", silent=TRUE)
    fitsim <- nlminb(objsim$par,objsim$fn,objsim$gr)
    as.list(sdreport(objsim),"Est")$gamma
}
sim <- replicate(100, doone())
dat <- odat
matplot(sim, type="l")
lines(as.list(sdreport(obj),"Est")$gamma, lwd=5, col="darkgreen")


