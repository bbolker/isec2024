library(scam)
library(glmmTMB)
library(marginaleffects)
library(ggplot2); theme_set(theme_bw())
library(purrr)
library(dplyr)
## what for tidy predictions? ggpredict, marginaleffects, emmeans?


## 1. simulate data from Holling type-2, type-3:
set.seed(101)

data("ReedfrogFuncresp", package = "emdbook")
dd <- ReedfrogFuncresp
## apropos("smooth.construct")
## ?smooth.construct.miso.smooth.spec
## scam smooth codes:
##  m = monotonic
##  p = p-spline
##  i/d = increasing/decreasing
##  cv = concavity

## linear fits
s1 <- scam(Killed/Initial ~ s(Initial, bs = "mpd"), data = dd)
s2 <- gam(Killed/Initial ~ s(Initial, bs = "tp", k = 8), data = dd)

## can't trust scam with binomial (yet)
## s3 <- scam(cbind(Killed, Initial-Killed) ~ s(Initial, bs = "mpd"), data = dd, family = binomial)
s3 <- glmmTMB(Killed/Initial ~ s(Initial, bs = "mpd"),
              weights = Initial, data = dd, family = binomial)

predfun <- function(m) {
    nd <- data.frame(Initial = 1:100)
    if (!inherits(m, "glmmTMB")) {
        c(marginaleffects::predictions(m, newdata = nd) |> as_tibble()
    } else {
        ## ugh ... get conditional predictions
        est <- predict(m, newdata = nd)
        ## X is 117 x 10, s$re.trans.U is 9x9 (why mismatch?)
        ## X has intercept??
    }
}

preds <- (list(scam_mpd = s1, gam_tp = s2, glmmTMB_mpd_binom = s3)
    |> map_dfr(predfun, .id = "model")
    |> select(model, Initial, prob = estimate, lwr = conf.low, upr = conf.high)
)
    

ggplot(preds, aes(Initial, prob)) +
    geom_line(aes(colour = model)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = model), colour = NA, alpha = 0.5) +
    expand_limits(y=0) +
    geom_point(data=dd, aes(y = Killed/Initial, size = Killed))

plot(Killed/Initial ~ Initial, data = dd, xlim = c(0, 100), ylim = c(0, 0.7))
with(pframe, lines(Initial, prob))
## confidence intervals?

plot(Killed ~ Initial, data = dd, xlim = c(0, 100))
lines(dd$Initial, predict(s1)*dd$Initial)
par(op)


## identity link
## non-integer # success warnings??  check these out.
f1 <- scam(cbind(Killed, Initial-Killed) ~ s(Initial, bs = "mpd"),
     family = binomial(),  ## try link == identity?
     data = dd)
plot(f1)
plot(predict(f1))

f2 <- scam(cbind(Killed, Initial-Killed) ~ s(Initial, bs = "mpd"),
     family = binomial(link = "identity"),  ## try link == identity?
     data = dd)
plot(predict(f2))

## no non-int success warnings ...
glm(cbind(Killed, Initial-Killed) ~ Initial,
     family = binomial(),  ## try link == identity?
     data = dd)

set.seed(7)
n <- 100
x <- c(-1,runif(n-1)*4-1) ## starting at -1 for a function to be zero at a start
z <- runif(n)
y <- exp(4*x)/(1+exp(4*x)) -0.01798621+ z*(1-z)*5 + rnorm(100)*.4
m1 <- scam(y~s(x,bs='miso')+s(z)) 
plot(m1,pages=1)
ewd <- data.frame(x=-1,z=0)
predict(m1,newd, type='terms')
