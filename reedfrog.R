library(scam)
library(glmmTMB)
library(marginaleffects)
library(ggplot2); theme_set(theme_bw())
library(purrr)
library(dplyr)
library(bbmle)
library(RTMB)
## what for tidy predictions? ggpredict, marginaleffects, emmeans?
source("funs.R")

## GOAL: fit Reed frog predation data, with CIs, using Holling type 2,
##  unrestricted GAM, shape-constrained GAM (both scam and RTMB), ...
##

## tests? parametric bootstrap, AIC, etc. ?
## (edf for RTMB fits?)

set.seed(101)

data("ReedfrogFuncresp", package = "emdbook")
dd <- ReedfrogFuncresp
ddx <- expand_bern(dd, response = "Killed", size = "Initial")

ddp0 <- data.frame(Initial = 1:100)
## works for scam/gam fits
pfun <- function(m) (c(marginaleffects::predictions(m, newdata = ddp0))
    |> as.data.frame()
    |> as_tibble()
    |> dplyr::select(prob = estimate, lwr = conf.low, upr = conf.high)
    |> mutate(Initial = ddp0$Initial, .before = 1)
)

## unconstrained GAM fit
m_gam_tp <- gam(cbind(Killed, Initial -Killed) ~ s(Initial, bs = "tp",
                                                          k = 8), data = dd,
                       family = binomial)
preds_gam_tp <- pfun(m_gam_tp)

## This doesn't work -- don't know if it's a thinko or just numerical nastiness ...
## holling  a*x/(1+a*h*x)
## initial slope = a; asymptote = 1/h ~ 0.5, 100
## holling prob = a/(1+a*h*x)
## holling via GLM; partial fractions!
## 1/prob = A + B/x = (Ax + B)/x ->
## prob = x/(Ax + B) = (1/A)/(1 + (B/A)*(1/x))
## a = 1/A; h = B

## 1/prob ~ b0 + b1/x
## prob = 1/(b0 + b1/x) = (1/b0)/(1/(1/b0) + (1/b0)*b1
if (FALSE) {
    ## don't want this lying around messing things up
    m_glm_holling <- glm(cbind(Killed, Initial -Killed) ~ I(1/Initial), data = dd,
                         family = binomial(link = "inverse"),
                         start = c(2, 0.02))
}
## coef(glm_holling)  ## bogus ...
## negative values for b1 ???
## plot(predict(glm_holling, newdata = data.frame(Initial = 1:100), type = "response"))


## fit Holling type 2 with mle2
m_mle2_holling <- bbmle::mle2(Killed ~ dbinom(prob = exp(loga)/(1+exp(loga)*exp(logh)*Initial),
                                                     size = Initial),
                                     start = list(loga = log(0.5), logh = log(0.01)),
                                     data = dd)

## predict number killed, divide by value
p0 <- predict(m_mle2_holling, newdata = list(Initial = 1:100))/(1:100)
rpars <- MASS::mvrnorm(1000, mu = coef(m_mle2_holling), Sigma = vcov(m_mle2_holling))
n_initial <- 1:100
preds <- apply(rpars, 1,
               function(x) { a <- exp(x[1]); h <- exp(x[2]); a/(1+a*h*n_initial)})
matplot(preds, type = "l", col = adjustcolor("black", alpha = 0.3), lty = 1)
ci <- t(apply(preds, 1, quantile, c(0.025, 0.975)))


preds_mle2_holling <- data.frame(Initial = 1:100, prob = p0, lwr = ci[,1], upr = ci[,2])

## fit Holling type 2 with RTMB
m_RTMB_holling <- fit_RTMB_holling2(dd)
preds_RTMB_holling <- predict_RTMB_holling2(dd, data.frame(Initial = 1:100),
                                         m_RTMB_holling$fit$par) |>
    dplyr::select(Initial, prob, lwr, upr)

## fit SCAM with RTMB (Laplace approx doesn't work, need random = NULL)
m_RTMB_mpd <- fit_mpd_fun(data = dd, response = "Killed",
                               size = dd$Initial, xvar = "Initial", family = "binomial", random = NULL)


m_RTMB_mpd$fit
ddp <- as.list(dd)
## don't go outside original range, technical issues with outer.ok in splineDesign ...
ddp$Initial <- c(dd$Initial, 5:100)
ddp$Killed <- c(dd$Killed, rep(NA_integer_, 96))
ddp <- as.data.frame(ddp) ## MUST have nrow()
k <- data.frame(Initial = smoothCon(s(Initial, bs="mpd"), data = dd, absorb.cons = TRUE)[[1]]$knots)
preds0 <- fit_mpd_fun(data = ddp, response = "Killed",
            size = ddp$Initial, xvar = "Initial", family = "binomial", random = NULL,
            knots = k, predict = TRUE,
            parms = with(m_RTMB_mpd$obj$env, parList(last.par.best)))
qq <- qnorm(0.975)
preds_RTMB_mpd <- (preds0
    |> filter(nm == "eta")
    |> slice_tail(n = 96)
    |> transmute(Initial = 5:100, prob = plogis(value), lwr = plogis(value-qq*sd),
                 upr = plogis(value+qq*sd))
)

m_scam_mpd <- scam(Killed ~ s(Initial, bs = "mpd"), data = ddx, family = binomial)

preds_scam_mpd <- pfun(m_scam_mpd)
(all_models <- ls(pattern="^m_"))
(all_preds <- ls(pattern="^preds_"))

## NOT 'preds_' (don't pollut
pred_frame <- (mget(all_preds)
    |> setNames(all_preds)
    |> bind_rows(.id = "model")
    |> mutate(across(model, ~gsub("preds_", "", .)))
)

ggplot(pred_frame, aes(Initial, prob)) +
    geom_line(aes(colour = model)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = model), colour = NA, alpha = 0.5) +
    expand_limits(y=0) +
    geom_point(data=dd, aes(y = Killed/Initial, size = Killed), alpha = 0.5) +
    facet_wrap(~model)
## CIs are not monotonic??


