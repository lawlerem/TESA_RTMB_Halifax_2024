# observations 
x <- rpois(1000,3)
ppois.u <- function(x, lambda){
  runif(length(x), ppois(x-1,lambda), ppois(x,lambda)) #uses the fact that ppois(-1,lambda)=0
}
U <- ppois.u(x,3)
Z <- qnorm(U)

