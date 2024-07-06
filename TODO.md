## to do

### reed frogs

* RTMB_mpd: penalization too strong? Is this (RE)ML vs GCV?
 * compare `m_scam_mpd[c("trA", "aic", "sp", "edf")]`
* get lambda, ecdf from RTMB_mpd, scam_mpd
* compare AIC values?
* ?? what does holling type 2 look like ??


## waterbugs

* fit!
* scam, RTMB (semimech)

## simulations?

## other

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

* Quinn and Deriso Table 3.2/p. 107 (search "Ricker" or "Table 3.2"):
https://books.google.ca/books?id=5FVBj8jnh6sC&printsec=frontcover&dq=quinn+deriso+quantitative+fisheries&hl=en&newbks=1&newbks_redir=0&sa=X&redir_esc=y#v=onepage&q=Ricker&f=false

* list of options?
* test 

scam smooth codes:
*  m = monotonic
*  p = p-spline
*  i/d = increasing/decreasing
*  cv = concavity

unimodal splines? (uniReg package)
