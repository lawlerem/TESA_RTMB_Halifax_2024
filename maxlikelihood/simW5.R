set.seed(123)
age<-sample(1:8, 123, replace=TRUE)
sex<-as.factor(sample(c("M","F"), 123, replace=TRUE))
k<-c(.4,.3)
Linf<-c(32,42)
logL<-log(Linf[sex]) + log(1-exp(-k[sex]*age))+rnorm(length(age), sd=0.05)
length<-exp(logL)
plot(age, length, col=sex, pch=4, cex=2, lwd=2)
legend("bottomright", legend=c("M", "F"), col=c("black","red"), pch=4, cex=2, lwd=2, lty=c(NA,NA), bty="n")

dat<-data.frame(age=age, sex=sex, length=round(length,1))

write.table(dat, file="files/length2.tab", row.names=FALSE, quote=FALSE)

