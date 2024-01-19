library(RTMB)
load("spatial_data.RData")

data<- rbind(
  data
)

correlation_function<- function(s1, s2, rho) {
  d<- sqrt( sum( (s1 - s2)^2 ) )
  C<- exp( -d/rho )
  return(C)
}
nll<- function(pars) {
  getAll(pars, data)
  sigmax<- exp(log_sigmax);  ADREPORT(sigmax)
  rho<- exp(log_rho); ADREPORT(rho)

  correlation_matrix<- matrix(0, nrow = length(x), ncol = length(x))
  for( i in seq_along(x) ) {
    for( j in seq_along(x) ) {
      correlation_matrix[i, j]<- correlation_function(
        data[i, c("lon", "lat")],
        data[j, c("lon", "lat")],
        rho
      )
      if( i != j ) correlation_matrix[i, j]<- 0.999 * correlation_matrix[i, j]
    }
  }

  x %~% dmvnorm(0, sigmax^2 * correlation_matrix)

  sigmay<- exp(log_sigmay)
  value<- OBS(value)
  value %~% dnorm(mu + x, sigmay)
}
pars<- list(
  log_sigmax = 0,
  log_rho = 0,
  mu = 0,
  x = numeric(nrow(data)),
  log_sigmay = 0
)
obj<- MakeADFun(
  nll,
  pars,
  random = "x"
)
opt<- nlminb(obj$par, obj$fn, obj$gr)
sdr<- sdreport(obj, opt$par)

obj$gr(opt$par)
# [1]  6.730522e-06 -2.776371e-07 -6.948028e-07 -1.288085e-06

set.seed(58236)
cc<- checkConsistency(obj, n = 500)
summary(cc)
# $joint
# $joint$p.value
# [1] 0.3798033

# $joint$bias
#    log_sigmax       log_rho            mu    log_sigmay
#  0.0045874021 -0.0013743910 -0.0045708001 -0.0001675774


# $marginal
# $marginal$p.value
# [1] 0.9605528

# $marginal$bias
#    log_sigmax       log_rho            mu    log_sigmay
#  0.0028865382 -0.0066614281 -0.0040109730  0.0001653636


logrho_profile<- TMB::tmbprofile(obj, "log_rho")
pdf("spatial_logrho_profile.pdf", width = 6, height = 4)
plot(logrho_profile, type = "l")
dev.off()

library(tmbstan)
fitmcmc2<- tmbstan(
  obj,
  chains = 1,
  iter = 10000,
  init = list(opt$par),
  lower = opt$par - 5 * summary(sdr)[1:length(opt$par), "Std. Error"],
  upper = opt$par + 20 * summary(sdr)[1:length(opt$par), "Std. Error"],
  laplace = TRUE
)

mc<- extract(
  fitmcmc2,
  pars = names(obj$par),
  inc_warmup = TRUE,
  permuted = FALSE
)
npar<- dim(mc)[3]

pdf("spatial_mcmc_chain.pdf", width = 10, height = 6)
layout(
  rbind(rep(1, npar),
  c(2:(npar + 1) ))
)
matplot(mc[, 1, ], type = "l")
for( i in 1:npar ) {
  pn<- attr(mc, "dimnames")$parameters[i]
  plot(density(mc[, 1, i]), main = pn, col = i, lwd = 3)
  abline(v = summary(sdr)[i, 1], lwd = 3)
  xx<- pretty(range(mc[, 1, i]), 100)
  lines(
    xx,
    dnorm(xx, summary(sdr)[i, 1], summary(sdr)[i, 2]),
    col = gray(.5),
    lwd = 3
  )
}
dev.off()

linear<- lm(mc[, , 2] ~ mc[, , 1])
pdf("spatial_mcmc_scatterplot.pdf", width = 6, height = 6)
plot(
  x = mc[, , 1],
  y = mc[, , 2],
  xlab = dimnames(mc)$parameters[[1]],
  ylab = dimnames(mc)$parameters[[2]],
  col = rgb(0, 0, 0, 0.1),
  pch = 19
)
abline(a = linear$coefficients[[1]], b = linear$coefficients[[2]], col = "red")
lines(
  c(opt$par["log_sigmax"], opt$par["log_sigmax"], -10),
  c(-10, opt$par["log_rho"], opt$par["log_rho"]),
  col = "red",
  lty = 2
)
dev.off()
