set.seed(123)
library(RTMB)
library(Matrix)

poly <- list(x = c(0.4, 0.37, 0.35, 0.34, 0.45, 0.47, 0.53, 0.64, 0.75,
                   0.79, 0.72, 0.75, 0.78, 0.81, 0.87, 0.78, 0.87, 0.74, 0.57, 0.48,
                   0.47, 0.28, 0.22, 0.28, 0.39, 0.4),
             y = c(0.48, 0.37, 0.16, 0.03, 0.08, 0.02, -0.04, 0.05, 0.01, 0.07, 0.33, 0.42, 0.41, 0.35,
                   0.49, 0.53, 0.65, 0.7, 0.67, 0.65, 0.61, 0.61, 0.44, 0.35, 0.46, 0.48))

plot(poly$x, poly$y, type="n")
polygon(poly, border=NA, col="lightblue")

DX<-.03
xgr<-seq(min(poly$x)-1.0e-9, max(poly$x)+DX, by=DX)
xcen<-xgr[-1]-0.5*DX

DY<-.03
ygr<-seq(min(poly$y)-1.0e-9, max(poly$y)+DY, by=DY)
ycen<-ygr[-1]-0.5*DY

CT<-matrix(rep(NA,((length(xgr)-1)*(length(ygr)-1))), ncol=(length(xgr)-1))

grd<-expand.grid(idx=1:ncol(CT), idy=1:nrow(CT))
idc<-which(sp::point.in.polygon(xcen[grd$idx],ycen[grd$idy],poly$x, poly$y)==1)
grd<-grd[idc,]

CT[cbind(grd$idy,grd$idx)]<-1:nrow(grd)

mycol <- rgb(0, 255, 0, max=255)

abline(v=xgr, col="gray")
abline(h=ygr, col="gray")
image(xcen, ycen, t(CT), add=T, col="cadetblue1")

xyc<-cbind(xcen[grd$idx],ycen[grd$idy], 1:nrow(grd))
text(xyc[,1], xyc[,2], labels=xyc[,3], cex=1.3, col="darkblue", font=2)

## find neighbors
AT<-array(c(CT,cbind(CT[,-1],NA), cbind(NA,CT[,-ncol(CT)]), rbind(CT[-1,],NA), rbind(NA,CT[-nrow(CT),])),dim=c(nrow(CT),ncol(CT),5))
nextTo<-apply(AT,3,function(x)x)
nextTo<-nextTo[!is.na(nextTo[,1]),]
nextTo<-nextTo[order(nextTo[,1]),]
D<-rowSums(!is.na(nextTo[,-1]), na.rm=TRUE)

Q0<-spMatrix(length(D), length(D))
diag(Q0)<-D
dummy<-apply(nextTo, 1, function(x){xx<-na.omit(x[-1]); cbind(rep(x[1],length(xx)),xx)})
nn<-do.call(rbind,dummy)
Q0[nn] <- -1

I<-spMatrix(nrow(Q0), nrow(Q0))
diag(I)<-1

# kind of data
dat<-list()
dat$x<-c(0.743, 0.725, 0.743, 0.721, 0.696, 0.725, 0.709, 0.74, 0.752,
0.801, 0.706, 0.705, 0.743, 0.677, 0.642, 0.679, 0.611, 0.732,
0.743, 0.57, 0.679, 0.586, 0.646, 0.594, 0.569, 0.501, 0.578,
0.488, 0.612, 0.428, 0.393, 0.381, 0.281, 0.306, 0.282, 0.368)

dat$y<-c(0.594, 0.554, 0.508, 0.172, 0.202, 0.536, 0.115, 0.527, 0.488,
0.459, 0.622, 0.505, 0.089, 0.116, 0.177, 0.303, 0.44, 0.601,
0.498, 0.234, 0.132, 0.404, 0.534, 0.569, 0.13, 0.186, 0.304,
0.531, 0.344, 0.331, 0.133, 0.528, 0.479, 0.584, 0.403, 0.091)

dat$obs<-c(83, 65, 64, 68, 80, 52, 59, 74, 50, 67, 44, 53, 44, 40, 41,
28, 38, 44, 36, 38, 29, 33, 31, 32, 21, 31, 22, 12, 10, 22, 16,
9, 9, 3, 3, 3)

xcut<-cut(dat$x,breaks=xgr)
ycut<-cut(dat$y,breaks=ygr)
dat$cellidx<-CT[cbind(as.integer(ycut), as.integer(xcut))]-1
dat$Q0<-Q0
dat$I<-I

par<-list()
par$logMu <- 0
par$logDelta <- 0
par$logSigma <- 0
par$logLambda <- rep(0, nrow(dat$Q0))

jnll <- function(par){
  getAll(par, dat)
  obs<-OBS(obs)
  delta <- exp(logDelta)
  sigma <- exp(logSigma)
  Q <- Q0 + delta*I
  ret <- 0
  ret <- -dgmrf(logLambda-logMu, Q=Q, scale=sigma, log = TRUE)
  ret <- ret-sum(dpois(obs, exp(logLambda[cellidx]), log=TRUE))
  out <- exp(logLambda)
  ADREPORT(out)
  ret
}

obj<-MakeADFun(jnll, par, random="logLambda")

fit <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj)
pl <- as.list(sdr,"Est")
plsd <- as.list(sdr,"Std")

print(as.list(sdr,what="Est", report=T)$out[1:5])

#stop()

zr<-c(0,ceiling(max(exp(pl$logLambda))))
cc<-colorRampPalette(c("blue","red"))(64)
lambda<-exp(pl$logLambda)
fitmat<-t(matrix(lambda[CT], nrow=nrow(CT), ncol=ncol(CT)))
fields::image.plot(xcen, ycen, fitmat, col=cc, xlab="", ylab="")
points(dat$x, dat$y, cex=dat$obs/max(dat$obs)*10, lwd=3)




library(tmbstan)
fitmcmc2<- tmbstan(
  obj,
  chains = 1,
  iter = 10000,
  init = list(fit$par),
  lower = fit$par - 20 * summary(sdr)[1:length(fit$par), "Std. Error"],
  upper = fit$par + 20 * summary(sdr)[1:length(fit$par), "Std. Error"],
  laplace = TRUE
)

mc<- extract(
  fitmcmc2,
  pars = names(obj$par),
  inc_warmup = TRUE,
  permuted = FALSE
)
npar<- dim(mc)[3]

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