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

## 1. simulate data from Holling type-2, type-3:
set.seed(101)

data("ReedfrogFuncresp", package = "emdbook")
dd <- ReedfrogFuncresp
ddx <- expand_bern(dd, response = "Killed", size = "Initial")
## apropos("smooth.construct")
## ?smooth.construct.miso.smooth.spec
## scam smooth codes:
##  m = monotonic
##  p = p-spline
##  i/d = increasing/decreasing
##  cv = concavity
## te = tensor
## d = double
## de = 'decreasing' (why not md?)

## unconstrained GAM fit
gam_rf_tp_binom <- gam(cbind(Killed, Initial -Killed) ~ s(Initial, bs = "tp",
                                                          k = 8), data = dd,

                       family = binomial)

## holling  a*x/(1+a*h*x)
## initial slope = a; asymptote = 1/h ~ 0.5, 100
## holling prob = a/(1+a*h*x)
## holling via GLM; partial fractions!
## 1/prob = A + B/x = (Ax + B)/x ->
## prob = x/(Ax + B) = (1/A)/(1 + (B/A)*(1/x))
## a = 1/A; h = B

## 1/prob ~ b0 + b1/x
## prob = 1/(b0 + b1/x) = (1/b0)/(1/(1/b0) + (1/b0)*b1

glm_rf_holling_binom <- glm(cbind(Killed, Initial -Killed) ~ I(1/Initial), data = dd,
                            family = binomial(link = "inverse"),
                            start = c(2, 0.02))
## coef(glm_rf_holling_binom)  ## bogus ...
## negative values for b1 ???
## plot(predict(glm_rf_holling_binom, newdata = data.frame(Initial = 1:100), type = "response"))


## fit Holling type 2 with mle2
mle2_rf_holling_binom <- bbmle::mle2(Killed ~ dbinom(prob = exp(loga)/(1+exp(loga)*exp(logh)*Initial),
                                                     size = Initial),
                                     start = list(loga = log(0.5), logh = log(0.01)),
                                     data = dd)
exp(coef(mle2_rf_holling_binom))
vcov(mle2_rf_holling_binom)

## fit Holling type 2 with RTMB
RTMB_rf_holling_binom <- fit_RTMB_holling2(dd)
pred <- pred_RTMB_holling2(dd, data.frame(Initial = 1:100), RTMB_rf_holling_binom$fit$par)

## fit SCAM with RTMB (Laplace approx doesn't work, need random = NULL)
RTMB_rf_scam_binom <- fit_mpd_fun(data = dd, response = "Killed", size = dd$Initial, xvar = "Initial", family = "binomial", random = NULL)

fit_mpd_fun(data = dd, response = "Killed", size = dd$Initial, xvar = "Initial", family = "binomial", random = NULL,
            smoothdata = dd, predict = TRUE,
            parms = with(RTMB_rf_scam_binom$obj$env, parList(last.par.best)))

scam_rf_mpd_binom <- scam(Killed ~ s(Initial, bs = "mpd"), data = ddx, family = binomial)

predictions(scam_rf_mpd_binom, newdata = data.frame(Initial = 1:100)) |> as_tibble()

predfun <- function(m) {
    nd <- data.frame(Initial = 1:100)
    if (!inherits(m, "glmmTMB")) {
        c(marginaleffects::predictions(m, newdata = nd) |> as_tibble())
    } else stop("can't do glmmTMB preds yet")
}

models <- ls(pattern="^(scam|gam|glm)_rf")
mod_list <- mget(models) |> setNames(models)
preds <- (mod_list
    |> map_dfr(predfun, .id = "model")
    |> select(model, Initial, prob = estimate, lwr = conf.low, upr = conf.high)
)
    

ggplot(preds, aes(Initial, prob)) +
    geom_line(aes(colour = model)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = model), colour = NA, alpha = 0.5) +
    expand_limits(y=0) +
    geom_point(data=dd, aes(y = Killed/Initial, size = Killed), alpha = 0.5) +
    facet_wrap(~model)
## CIs are not monotonic??


scam_pos <- match("package:scam", search())
aa <- apropos("smooth.construct", where = TRUE)
scam_smooths <- unname(aa[names(aa) == scam_pos]) |>
    gsub(pattern = "[.]?smooth\\.(construct|spec)[.]?", replacement = "")
