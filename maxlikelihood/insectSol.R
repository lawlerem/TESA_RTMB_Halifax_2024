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
obj <- MakeADFun(f, par, silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
 
map1=list(logAlpha=factor(c(1,1,2,3,4,1)))
obj1 <- MakeADFun(f, par, map=map1, silent=TRUE)
opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)
1-pchisq(2*(opt1$obj-opt$obj),2)
#
#  0.3982677
#
map2=list(logAlpha=factor(c(NA,NA,2,3,4,NA)))
par2<-par
par2$logAlpha[c(1,2,6)]<-log(15)
obj2 <- MakeADFun(f, par2, map=map2, silent=TRUE)
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
1-pchisq(2*(opt2$obj-opt1$obj),1)
#
# 0.4410911
#
