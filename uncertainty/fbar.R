th <- c(-1.13, -.75, -.94)
V <- matrix(c(0.0222, 0.0135, 0.0114, 0.0135, 0.0169, 0.0137, 0.0114, 0.0137, 0.0191), 3)


# by hand 
grad <- exp(th)/sum(exp(th))
logFbar <- log(mean(exp(th)))
sdlogFbar <- sqrt(grad%*%V%*%grad)[1,1]
# sim validation 
# sd(apply(MASS::mvrnorm(1000, th, V), 1, function(x) log(mean(exp(x)))))
cilogFbar <- logFbar+c(-2,2)*sdlogFbar
ciFbar1 <- exp(cilogFbar)

library(RTMB)
mymean <- function(x)sum(x)/length(x)           ### annoying talk to kkr 
calcFbar <- function(th)log(mymean(exp(th)))   
F <- MakeTape(calcFbar,th)
est <- F(th)
grad <- F$jacobian(th)
sd <- sqrt(grad%*%V%*%t(grad))[1,1]
cilogFbar <- est+c(-2,2)*sd
ciFbar2 <- exp(cilogFbar)

