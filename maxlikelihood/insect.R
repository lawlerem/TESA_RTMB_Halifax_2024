library(RTMB)

# for data we use the built-in data "InsectSprays"
par <- list(logAlpha=rep(0,nlevels(InsectSprays$spray)))
f<-function(par){
  getAll(InsectSprays, par)
  nll <- 0
  for(i in 1:length(count)){
    lambda <- exp(logAlpha[spray[i]])
    nll <- nll - dpois(count[i],lambda,log=TRUE)
  }
  nll
}
obj <- MakeADFun(f, par)
opt <- nlminb(obj$par, obj$fn, obj$gr)
