library(RTMB)
load("age_length_data.RData")
load("age_length_helpers.RData")
source("age_length.R")

mesh_info<- list(
  mesh = mesh,
  spde = spde,
  interpolator_data = interpolator_data,
  interpolator_prediction = interpolator_prediction
)

par<- list(
  log_L_inf = log(max(cod$length)),
  log_tau = 0,
  log_kappa = 0,
  log_rate = numeric(mesh$n),
  mu = 0,
  nlog_t0 = 0,
  log_obs_sd = 0
)

nll<- function(par) {
  getAll(par, cod, mesh_info)
  L_inf<- exp(log_L_inf)

  # GMRF prior for rate
  tau<- exp(log_tau)
  kappa<- exp(log_kappa)
  Q<- tau * (kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2)
  log_rate %~% dgmrf(mu, Q)

  rate<- exp(interpolator_data %*% log_rate)
  ADREPORT(rate)
  t0<- -exp(nlog_t0)
  pred<- log(L_inf) + log(1 - exp(-rate * (age - t0)))
  ADREPORT(pred)

  predictions<- exp(interpolator_prediction %*% log_rate)
  ADREPORT(predictions)

  obs_sd<- exp(log_obs_sd)
  log(length) %~% dnorm(
    pred,
    obs_sd
  )
}
obj<- MakeADFun(nll, par, random = "log_rate")
spatial_fit<- nlminb(obj$par, obj$fn, obj$gr)
sdrep<- sdreport(obj)
pl<- as.list(sdrep, "Est", report = TRUE)
plsd<- as.list(sdrep, "Std", report = TRUE)


prediction_locations<- cbind(
  rate_est = pl$predictions,
  rate_se = plsd$predictions,
  prediction_locations
)
prediction_locations[rowSums(interpolator_prediction) == 0, c("rate_est", "rate_se")]<- NA

prediction_locations$lower<- prediction_locations$rate_est - 1.96 * prediction_locations$rate_se
prediction_locations$upper<- prediction_locations$rate_est + 1.96 * prediction_locations$rate_se



color_range<- seq(
  min(prediction_locations$rate_est, na.rm = TRUE),
  max(prediction_locations$rate_est, na.rm = TRUE),
  length.out = 10
)

pdf("age_length_predictions.pdf", width = 4, height = 4)
plot(
  x = prediction_locations[, "lon"],
  y = prediction_locations[, "lat"],
  col = "black",
  bg = hcl.colors(10, "Blues", rev = TRUE)[
    cut(
      prediction_locations$rate_est,
      color_range,
      include.lowest = TRUE
    )
  ],
  ,
  xlab = "lon",
  ylab = "lat",
  main = "Mean",
  pch = 21
)
dev.off()




color_range<- seq(
  min(prediction_locations$lower, na.rm = TRUE),
  max(prediction_locations$upper, na.rm = TRUE),
  length.out = 10
)

# Darker areas = higher values
pdf("age_length_predictions_CI.pdf", width = 8, height = 6)
par(mfrow = c(1, 3))
plot(
  x = prediction_locations[, "lon"],
  y = prediction_locations[, "lat"],
  col = "black",
  bg = hcl.colors(10, "Blues", rev = TRUE)[
    cut(
      prediction_locations$lower,
      color_range,
      include.lowest = TRUE
    )
  ],
  xlab = "lon",
  ylab = "lat",
  main = "Lower 95%",
  pch = 21
)
plot(
  x = prediction_locations[, "lon"],
  y = prediction_locations[, "lat"],
  col = "black",
  bg = hcl.colors(10, "Blues", rev = TRUE)[
    cut(
      prediction_locations$rate_est,
      color_range,
      include.lowest = TRUE
    )
  ],
  ,
  xlab = "lon",
  ylab = "lat",
  main = "Mean",
  pch = 21
)
plot(
  x = prediction_locations[, "lon"],
  y = prediction_locations[, "lat"],
  col = "black",
  bg = hcl.colors(10, "Blues", rev = TRUE)[
    cut(
      prediction_locations$upper,
      color_range,
      include.lowest = TRUE
    )
  ],
  xlab = "lon",
  ylab = "lat",
  main = "Upper 95%",
  pch = 21
)
dev.off()


# Likelihood ratio test to determine if we can treat growth rate as spatially constant
spatial_nll<- spatial_fit$objective
non_spatial_nll<- non_spatial_fit$objective
df<- length(spatial_fit$par) - length(non_spatial_fit$par)
pchisq(-2 * (spatial_nll - non_spatial_nll), df = df, lower.tail = FALSE)
# p ~ 0.095 (what do we conclude?)
