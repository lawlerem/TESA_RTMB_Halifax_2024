x <- c(73, 75, 84, 76, 93, 79, 85, 80, 76, 78, 80)

o <- order(x)
sigma <- sd(x)

# a) all observed 

la<-function(mu) -sum(dnorm(x, mu, sd=sigma, log=TRUE))

plot(Vectorize(la), 60,100)

# b) only mean

lb<-function(mu) -dnorm(mean(x), mu, sd=sigma/sqrt(11), log=TRUE)

plot(lb, 60,100)

# c) only median

lc<-function(mu) -dnorm(median(x), mu, sd=sigma, log=TRUE)-5*pnorm(median(x),mu,sigma,log=TRUE)-5*pnorm(median(x),mu,sigma,log=TRUE, lower.tail=FALSE)

plot(lc, 60,100)

# d) only min and max

ld<-function(mu) -dnorm(min(x), mu, sd=sigma, log=TRUE)-dnorm(max(x), mu, sd=sigma, log=TRUE)-9*log(pnorm(max(x),mu,sigma)-pnorm(min(x),mu,sigma))

plot(ld, 60,100)

# e) two lowest

le<-function(mu) -dnorm(x[o[1]], mu, sd=sigma, log=TRUE)-dnorm(x[o[2]], mu, sd=sigma, log=TRUE)-9*pnorm(x[o[2]],mu,sigma, log=TRUE, lower.tail=FALSE)

plot(le, 60,100)

xx<-seq(60,100, length=1000)
l<-cbind(A=Vectorize(la)(xx), B=lb(xx), C=lc(xx), D=ld(xx), E=le(xx))
l<-apply(l,2,function(x)x-min(x))
matplot(xx, l, type="l", lwd=3, lty="solid", col=c("darkred", "darkgreen","darkblue","orange","hotpink"))
legend("top", legend=c("A","B","C","D","E"), col=c("darkred", "darkgreen","darkblue","orange","hotpink"), lwd=3, bty="n")
