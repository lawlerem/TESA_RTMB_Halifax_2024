library ( RTMB )
dat <- list (X=2)
par <- list ( alpha =0)

f <- function ( par ){
  p <- plogis ( par $ alpha )                         ### Notice exp ( alpha )/ (1+ exp ( alpha ))
  -dbinom ( dat $X ,100 ,p, log = TRUE )
}

obj <- MakeADFun (f, par , silent = TRUE )
opt <- nlminb ( obj $par , obj $fn , obj $gr)
sdr <- sdreport ( obj )
summary (sdr)
#       Estimate Std . Error
# alpha -3.89182   0.7142857
pl <-as. list (sdr ," Est ")
plsd <-as. list (sdr ," Std ")
plogis (pl$ alpha +c( -2 ,2)* plsd $ alpha )
# 0.004867034 0.078475060                             ### Use same transformation to calculate CI
