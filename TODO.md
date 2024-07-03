## to do

* simulations
* play with reed frog data from emdbook
* RTMB?? mixed models??
* dig out McCoy et al data
* simulation-based examples
* hypothesis tests (goodness-of-fit/AIC, monotonicity, concavity, etc.)
* dynamical sensitivity???

## technical bits

* plotting? `gratia`?
* figure out weights for scam: are weights precision weights or analytic weights (sensu Lumley https://notstatschat.rbind.io/2020/08/04/weights-in-statistics/) ?

  are weights invariant under scaling (precision weights) or not (frequency weights)?

## waterbugs

* what did we do originally?
* functions for attack(size): Ricker, power-Ricker, logistic, hyperbolic, exponential,
* functions for handling(size): exponential, linear, proportional, independent

prob = 1/((1/a) + h*init_dens)

* gam, scam, glmmTMB, RTMB, JAGS ?
*  x simulated data, reedfrog, McCoy data
*  x Gaussian, binomial
*  x GAM, scGAM, parametric

regularization??

## talk outline

* semipar models
* Levins
* shape-constrained models
* `gam`, `scam`
* bases and constraints
* test 

scam smooth codes:
*  m = monotonic
*  p = p-spline
*  i/d = increasing/decreasing
*  cv = concavity

unimodal splines? (uniReg package)
