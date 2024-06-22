library(emdbook)
library(scam)
data("ReedfrogFuncresp")
apropos("smooth.construct")
?smooth.construct.miso.smooth.spec
## scam smooth codes:
##  m = monotonic
##  p = p-spline
##  i/d = increasing/decreasing
##  cv = concavity
odd <- ReedfrogFuncresp
plot(Killed/Initial ~ Initial, data = dd, xlim = c(0, 100))
s1 <- scam(Killed/Initial ~ s(Initial, bs = "mpd"), data = dd)
lines(dd$Initial, predict(s1))

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
newd <- data.frame(x=-1,z=0)
predict(m1,newd, type='terms')
