## make sure we have up-to-date versions
remotes::install_github("bbolker/reformulas")
remotes::install_github("glmmTMB/glmmTMB/glmmTMB")

library(glmmTMB)
library(mgcv)
library(scam)
library(Matrix)
library(RTMB)
library(cowplot)

## simulated example
set.seed(101)
dd <- data.frame(x=seq(-5, 5, length = 101))
dd <- within(dd, {
             mu <- 1 + x - x^3/4
             y <- rnorm(length(mu), mu, sd = 1)
             })

par(las = 1, bty = "l") ## cosmetic
colvec <- c(2,4:6, 8:9)

m_gam_tp_gcv <- gam(y ~ s(x, bs = "tp"), method="GCV.Cp", data = dd)
m_gam_tp_reml <- gam(y ~ s(x, bs = "tp"), method = "REML", data = dd)

m_scam_mpd_gcv <- scam(y ~ s(x, bs = "mpd"), data = dd)
m_scam_tp_gcv <- scam(y ~ s(x, bs = "tp"), data = dd)

m_glmmTMB_mpd_reml <- glmmTMB(y ~ s(x, bs = "mpd"), data = dd, REML = TRUE)
m_glmmTMB_tp_reml <- glmmTMB(y ~ s(x, bs = "tp"), data = dd, REML = TRUE)

nm0 <- ls(pattern = "m_.*_.*_.*")
nm <- nm0 |> gsub(pattern = "^m_", replacement = "")

predmat <- mget(nm0) |> lapply(predict) |> do.call(what=cbind)
colnames(predmat) <- nm

focal_col <- "scam_mpd_gcv"

predmat_diff <- sweep(predmat[,nm != focal_col], 1, FUN = "-", predmat[, focal_col])
matplot(predmat_diff, col = colvec, type = "l")
legend("topleft", legend = setdiff(nm, focal_col), col = colvec, lty = 1:8)
## gam_tp_gcv = scam_tp_gcv
## gam_tp_reml = glmmTMB_tp_reml
## glmmTMB_mpd_reml is just bogus

## note, REML appears to be much closer to GCV for 'ps' basis (instead of 'tp')

######
## glmmTMB ignores need to exponentiate parameters. What else?

sm2 <- smoothCon(s(x, bs = "ps"), data = dd, absorb.cons = TRUE)[[1]]
sm0 <- smoothCon(s(x, bs = "tp"), data = dd, absorb.cons = TRUE)[[1]]
sm1 <- smoothCon(s(x, bs = "mpd"), data = dd, absorb.cons = TRUE)[[1]]

ifun <- function(x) image(Matrix(x), sub = "",
                          xlab = "", ylab = "")

p1 <- ~matplot(sm0$X, type = "l", main = "bs = 'tp'")
p2 <- ~matplot(sm1$X, type = "l", main = "bs = 'mpd'")
p3 <- ~matplot(sm2$X, type = "l", main = "bs = 'ps'")
p4 <- ifun(sm0$S[[1]])
p5 <- ifun(sm1$S[[1]])
p6 <- ifun(sm2$S[[1]])
plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2)
## not sure exactly why S is simpler for mpd than tp ... ?

## smooth2random messes things up by reparameterizing -- we need original
##  b values for deciding where to apply constraints/exponentiation

## what is in a tp (gam) smooth object vs the mpd object?
setdiff(names(sm0), names(sm1))
## drop.null, UZ, xU, shift
setdiff(names(sm1), names(sm0))
## extra stuff ...
## cmX: original column centers
## Sigma: cumulation matrix (D in Pya and Wood)
## P: square root of S
## p.ident: which parameters have positivity constraint?
## C: constraint matrix (empty)
## Xdf1, Xdf2: model matrices for 1st and 2d derivatives
## knots
## m: order
## [1] "cmX"     "Sigma"   "P"       "p.ident" "C"       "Xdf1"    "Xdf2"   
## [8] "knots"   "m"

## n.b. X already includes multiplication by D

## contents of smoothCon
## If constraints are
## to be absorbed then the objects will have attributes ‘"qrc"’ and
## ‘"nCons"’. ‘"nCons"’ is the number of constraints. ‘"qrc"’ is
## usually the qr decomposition of the constraint matrix (returned by
## qr), but if it is a single positive integer it is the index of
## the coefficient to set to zero, and if it is a negative number
## then this indicates that the parameters are to sum to zero.

## implement in RTMB

## make b1 parameters non-zero so we can see their effect when
## b1 is fixed parameter vec
set.seed(101)
parameters <- list(
    b0 = 0,
    ## b1 = rep(0, length(sm1$p.ident)),
    b1 = rnorm(length(sm1$p.ident)),
    log_smSD = 0,
    log_rSD = 0
)

## shouldn't need log(det(S)) as it's constant
tmbdat_mpd1 <- c(as.list(dd), list(p.ident = sm1$"p.ident", S = sm1$S[[1]], X = sm1$X))

## dmvnorm with rank-deficient covariance matrix

## not working yet, should try with tp() (and smooth2random)
## can't use %~% format if we want to add a penalty
f_mpd1 <- function(parms) {
    getAll(tmbdat, parms)
    b_pos <- exp(b1)
    mu <- b0 + X %*% b_pos
    nll <- -1*sum(dnorm(y, mu, exp(log_rSD), log = TRUE))
    ## lambda = 1/sigma_sm^2
    ## MVgauss NLL = (1/2) (n*log(2pi) + log(det(Sigma)) + bT Sigma^{-1} b)
    ## Sigma = sd^2 S^{-1}
    ## log(det(Sigma)) = 2*logsd - log(det(S))
    ## MVNLL = C + logsd + 1/sd^2 (bT S b)
    pen <- (exp(-2*log_smSD) * (t(b_pos) %*% S %*% b_pos) + 2*log_smSD)/2
    REPORT(mu)
    nll + pen
}

obj_mpd1 <- MakeADFun(f_mpd1, parameters,
                      random=c("b1"),
                      silent = TRUE)
res_mpd1 <- with(obj_mpd1, nlminb(par, fn, gr))
plot(dd$x, obj_mpd1$report()$mu, type = "l")
points(dd$x, dd$y)
lines(dd$x, predict(m_scam_mpd_gcv), col = 2)

####

tmbdat <- c(as.list(dd), list(p.ident = sm1$"p.ident", S = sm1$S[[1]], X = sm1$X))
f <- function(parms) {
    getAll(tmbdat, parms)
    mu <- b0 + X %*% b1
    nll <- -1*sum(dnorm(y, mu, exp(log_rSD), log = TRUE))
    pen <- (exp(-log_smSD) * (t(b1) %*% S %*% b1) + log_smSD)/2
    REPORT(mu)
    nll + pen
}

## dbeta.rho[, j] <- -sp[j] * P %*% (t(P) %*% (S[[j]] %*% 
##                  beta))

er <- eigen(sm1$S[[1]], symmetric=TRUE)

er$values <- zapsmall(er$values)
rS <- crossprod(sqrt(sqrt(er$values))*t(er$vectors))
b1 <- parameters$b1
stopifnot(all.equal(sum((rS %*% b1)^2), c(t(b1) %*% sm1$S[[1]] %*% b1)))

p2 <- parameters
p2$log_smSD <- p2$log_smSD+0.01
f(parameters)
f(p2)


obj$fn()
obj$gr()
res1 <- with(obj, optim(par, fn, method = "Nelder-Mead", control = list(maxit = 1000, trace = 100)))
res2 <- with(obj, nlminb(par, fn, gr,
             lower = c(-Inf, rep(0, 9), -Inf, -Inf), control = list(eval.max = 500)))

res2 <- with(obj, nlminb(par, fn, gr), control = list(eval.max = 500))
## looks like tp is fine ... ???
plot(tmbdat$x, obj$report()$mu, type = "l")
lines(dd$x, predict(m_scam_mpd_gcv))
points(tmbdat$x, tmbdat$y)
res1$par
res2$par
m_scam_mpd_gcv


plot(dd$x, predict(m_scam_mpd_gcv))
lines(dd$x, predict(m_scam_tp_gcv), col = 2, lwd = 2)
with(dd, points(x, y, col = adjustcolor("black", alpha = 0.3)))
lines(dd$x, obj$report()$mu)

## try to make it match with a fixed smoothing penalty???
