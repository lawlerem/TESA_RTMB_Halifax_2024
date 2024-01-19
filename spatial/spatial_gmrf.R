library(RTMB)
load("spatial_data.RData") # gets data, prediction_locations, raster
load("spatial_helpers.RData") # gets mesh, spde, interpolators

mesh_info<- list(
  mesh = mesh,
  spde = spde,
  interpolator_data = interpolator_data,
  interpolator_prediction = interpolator_prediction
)

nll <- function(par) {
  getAll(par, data, mesh_info)
  tau <- exp(log_tau)
  kappa <- exp(log_kappa)
  # Computation of precision matrix Q here is a standard line
  # of code that you can copy/paste into your own code
  Q <- tau * (kappa^4 * spde$c0 + 2 * kappa^2 * spde$g1 + spde$g2)
  x %~% dgmrf(0, Q)

  data_predictions <- mu + interpolator_data %*% x
  obs_sd <- exp(log_obs_sd)
  value %~% dnorm(data_predictions, obs_sd)

  predictions <- mu + interpolator_prediction %*% x
  ADREPORT(predictions)
}

par <- list(
  log_tau = 0,
  log_kappa = 0,
  x = numeric(mesh$n),
  mu = 0,
  log_obs_sd = 0
)
obj <- MakeADFun(nll, par, random = "x")
fit <- nlminb(obj$par, obj$fn, obj$gr)
sdrep <- sdreport(obj)
pl <- as.list(sdrep, "Estimate", report = TRUE)
plsd <- as.list(sdrep, "Std. Error", report = TRUE)

est_field <- matrix(pl$predictions, nrow = nrow(raster))
sd_field <- matrix(plsd$predictions, nrow = nrow(raster))

par(mfrow = c(3, 1))
image(raster)
image(est_field)
image(sd_field)
