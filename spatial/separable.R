library(RTMB)
load("matrix_ar1.RData")

nll<- function(pars) {
  getAll(pars)
  rho<- 2 * plogis(working_rho) - 1
  sigma<- exp(log_sigma)

  ar1_within_rows<- function(x) dautoreg(
    x,
    phi = rho[1],
    log = TRUE,
    scale = sigma[1]
  )
  ar1_within_cols<- function(x) dautoreg(
    x,
    phi = rho[2],
    log = TRUE,
    scale = sigma[2]
  )
  dmatrix_ar1<- dseparable(ar1_rows, ar1_cols)
  x %~% dmatrix_ar1()

  # Can also scale x directly by doing...
  # ar1_rows<- function(x) dautoreg(x, phi = rho[1], log = TRUE)
  # ar1_cols<- function(x) dautoreg(x, phi = rho[2], log = TRUE)
  # dmatrix_ar1<- dseparable(ar1_along_rows, ar1_along_columns)
  # x %~% dmatrix_ar1(scale = sigma[1] * sigma[2])

  # You can input a scaling matrix...
  # scale_matrix<- matrix(sigma[1] * sigma[2], nrow(x), ncol(x))
  # x %~% dmatrix_ar1(scale = scale_matrix)

  # Can write the negative likelihood directly with...
  # nll<- -dseparable(ar1_along_rows, ar1_along_columns)(x)

  y %~% dpois(exp(mu + x))
}
pars<- list(
  mu = 0,
  working_rho = qlogis(0.5 * c(0.1, 0.9) + 0.5),
  log_sigma = log(c(0.3, 5)),
  x = matrix(
    0,
    nrow(y),
    ncol(y)
  )
)
obj<- MakeADFun(nll, pars, random = "x")
opt<- nlminb(obj$par, obj$fn, obj$gr)
