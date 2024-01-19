library(RTMB)
load("spatial_data.RData")

data<- rbind(
  data,
  prediction_locations
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
      # Ensure the matrix is invertible
      if( i != j ) correlation_matrix[i, j]<- 0.999 * correlation_matrix[i, j]
    }
  }

  x %~% dmvnorm(mu, sigmax^2 * correlation_matrix)

  sigmay<- exp(log_sigmay)
  value<- OBS(value)
  value[!is.na(value)] %~% dnorm(x[!is.na(value)], sigmay)
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
sdr_est<- as.list(sdr, "Est")
sdr_se<- as.list(sdr, "Std")


sigmax2_CI<- exp(
  sdr_est$log_sigmax + c(-1.96, 1.96) * sdr_se$log_sigmax
)^2
rho_CI<- exp(
  sdr_est$log_rho + c(-1.96, 1.96) * sdr_se$log_rho
)


prediction_locations$value<- sdr_est$x[is.na(data$value)]
prediction_locations$se<- sdr_se$x[is.na(data$value)]

pdf("spatial_prediction.pdf", width = 8, height = 4)
par(mfrow = c(1, 2))
image(raster, main = "Predicted Value")
points(
  x = prediction_locations$lon,
  y = prediction_locations$lat,
  col = "black",
  bg = hcl.colors(12, "YlOrRd", rev = TRUE)[
    cut(
      prediction_locations$value,
      quantile(
        c(raster),
        seq(0, 1, length.out = 12)
      ),
      include.lowest = TRUE
    )
  ],
  cex = 1.2,
  pch = 21,
  xlim = c(0, 1),
  ylim = c(0, 1)
)
points(
  x = data$lon[!is.na(data$value)],
  y = data$lat[!is.na(data$value)],
  col = "black",
  cex = 0.3,
  pch = 19
)
image(raster, main = "Standard Error")
points(
  x = prediction_locations$lon,
  y = prediction_locations$lat,
  col = "black",
  bg = hcl.colors(5, "Blues", rev = TRUE)[
    cut(
      prediction_locations$se,
      seq(0, max(prediction_locations$se), length.out = 5),
      include.lowest = TRUE
    )
  ],
  cex = 1.2,
  pch = 21,
  xlim = c(0, 1),
  ylim = c(0, 1)
)
points(
  x = data$lon[!is.na(data$value)],
  y = data$lat[!is.na(data$value)],
  col = "black",
  cex = 0.3,
  pch = 19
)
dev.off()
