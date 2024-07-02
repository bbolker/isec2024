remotes::install_github("bbolker/reformulas")
remotes::install_github("glmmTMB/glmmTMB/glmmTMB")

library(glmmTMB)
library(mgcv)
library(scam)
library(Matrix)
library(RTMB)
set.seed(101)
dd <- data.frame(x=seq(-5, 5, length = 101))
dd <- within(dd, {
             mu <- 1 + x - x^3/4
             y <- rnorm(length(mu), mu, sd = 1)
             })

par(las = 1, bty = "l") ## cosmetic
colvec <- c(2,4:6, 8:9)

m_gam_tp_gcv <- gam(y ~ s(x, bs = "tp"), data = dd)
m_gam_tp_reml <- gam(y ~ s(x, bs = "tp"), method = "REML", data = dd)

m_scam_mpd_gcv <- scam(y ~ s(x, bs = "mpd"), data = dd)
m_scam_tp_gcv <- scam(y ~ s(x, bs = "tp"), data = dd)

m_glmmTMB_mpd_reml <- glmmTMB(y ~ s(x, bs = "mpd"), data = dd, REML = TRUE)
m_glmmTMB_tp_reml <- glmmTMB(y ~ s(x, bs = "tp"), data = dd, REML = TRUE)


predmat <- cbind(predict(m_gam_tp_gcv),
                 predict(m_gam_tp_reml),
                 predict(m_scam_mpd_gcv),
                 predict(m_scam_tp_gcv),
                 predict(m_glmmTMB_mpd_reml),
                 predict(m_glmmTMB_tp_reml))

predmat_diff <- sweep(predmat[,-3], 1, FUN = "-", predmat[,3])
matplot(predmat_diff)

sm0 <- smoothCon(s(x, bs = "tp"), data = dd, absorb.cons = TRUE)[[1]]
sm1 <- smoothCon(s(x, bs = "mpd"), data = dd, absorb.cons = TRUE)[[1]]

## smooth2random messes things up by reparameterizing -- we need original
##  b values for deciding where to apply constraints/exponentiation
smooth2random(sm1, vnames = "", type = 2)
names(sm1)
setdiff(names(sm0), names(sm1))
## drop.null, UZ, xU, shift
setdiff(names(sm1), names(sm0))
## [1] "cmX"     "Sigma"   "P"       "p.ident" "C"       "Xdf1"    "Xdf2"   
## [8] "knots"   "m"

## p.ident is an indicator of which coefficients must be positive (exponentiated)
## Xdf1, Xdf2 are model matrices for 1st and 2d derivatives
sm1$p.ident


## contents of smoothCon
## If constraints are
## to be absorbed then the objects will have attributes ‘"qrc"’ and
## ‘"nCons"’. ‘"nCons"’ is the number of constraints. ‘"qrc"’ is
## usually the qr decomposition of the constraint matrix (returned by
## qr), but if it is a single positive integer it is the index of
## the coefficient to set to zero, and if it is a negative number
## then this indicates that the parameters are to sum to zero.

## implement in RTMB

## rearrange S/b into Xf/Xr

## X does have D built into it (from Pya and Wood)
## debug(smooth.construct.mpd.smooth.spec)

parameters <- list(
    b0 = 0,
    b1 = rep(0, length(sm1$p.ident)),
    log_smSD = 0,
    log_rSD = 0
)

## shouldn't need log(det(S)) as it's constant
tmbdat <- c(as.list(dd), list(p.ident = sm$"p.ident", S = sm$S[[1]], X = sm$X))

## dmvnorm with rank-deficient covariance matrix

## not working yet, should try with tp() (and smooth2random)
f <- function(parms) {
    getAll(tmbdat, parms)
    y <- OBS(y)
    mu <- b0 + X %*% exp(b1)
    nll <- 0
    y %~% dnorm(mu, exp(log_rSD)))
    pen <- (exp(log_smSD) * (t(b1) %*% S %*% b1) + log_smSD)/2
    REPORT(mu)
    nll + pen
}

f(parameters)

obj <- MakeADFun(f, parameters,
                 random=c("b1"),
                 silent = TRUE)

obj$fn()
obj$gr()
res <- with(obj, nlminb(par, fn, gr), control = list(eval.max = 500))

m_scam_mpd_gcv


plot(dd$x, predict(m_scam_mpd_gcv))
lines(dd$x, predict(m_scam_tp_gcv), col = 2, lwd = 2)
with(dd, points(x, y, col = adjustcolor("black", alpha = 0.3)))
lines(dd$x, obj$report()$mu)
