library(RTMB)
load("age_length_data.RData")
cod$log_length<- log(cod$length)

par<- list(
  log_L_inf = log(max(cod$length)),
  log_rate = 0,
  nlog_t0 = 0,
  log_obs_sd = 0
)
nll<- function(par) {
  getAll(par, cod)
  L_inf<- exp(log_L_inf)
  rate<- exp(log_rate)
  t0<- -exp(nlog_t0)
  pred<- log(L_inf) + log(1 - exp(-rate * (age - t0)))
  ADREPORT(pred)

  obs_sd<- exp(log_obs_sd)
  ADREPORT(obs_sd)
  log_length<- OBS(log_length)
  log_length %~% dnorm(
    pred,
    obs_sd
  )
}
obj<- MakeADFun(nll, par)
non_spatial_fit<- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport (obj)
pl <- as.list (sdr, "Est")
plsd <- as.list (sdr, "Std")
plr <- as.list (sdr, "Est", report = TRUE)
plrsd <- as.list (sdr, "Std", report = TRUE)
plot(cod[, c("age", "length")], cex.axis = 1.7, cex.lab = 1.7)
o <-order(cod$age)
lines (cod$age[o], exp(plr$pred[o]), lwd = 2, col = "orange")
lines (cod$age[o], exp(plr$pred[o] - 2 * plrsd$pred[o]), lwd = 2, col = "orange", lty = 2)
lines (cod$age[o], exp(plr$pred[o] + 2 * plrsd$pred[o]), lwd = 2, col = "orange", lty = 2)
lines (cod$age[o], exp(plr$pred[o] - 2 * (plrsd$pred[o] + plr$obs_sd[[1]])), lwd = 2, col = "orange", lty = 3)
lines (cod$age[o], exp(plr$pred[o] + 2 * (plrsd$pred[o] + plr$obs_sd[[1]])), lwd = 2, col = "orange", lty = 3)
legend(
  "topleft",
  col = "orange",
  lty = c(1, 2, 3),
  legend = c(
    "Estimate",
    "95% CI for Mean",
    "95% CI for Individuals"
  ),
  cex = 1.7
)
