## to do

##  simulations

## reed frogs

* RTMB_mpd: penalization too strong? Is this (RE)ML vs GCV?
 * compare `m_scam_mpd[c("trA", "aic", "sp", "edf")]`
 * works OK for binom_test (better data!)
 
* get lambda, ecdf from RTMB_mpd, scam_mpd
   * `sp` is multiplied by 
* compare AIC values?
* ?? what does holling type 2 look like ??

* tmbstan
* constrained optimization?
* bad knots? too many, wrong place? (hard to adjust ...)

## waterbugs

* fit!
* pix
* AIC tables for various parametric fits
* scam, RTMB (semimech)?? one smooth for mpd, 


larvae of red-eyed treefrogs Agalychnis callidryas and two species of aquatic invertebrate predators, adult predatory water bugs (Belostoma sp. Belostomatidae) and dragonfly nymphs (Pantala flavescens Libellulidae)


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

* list of options?
* test 

scam smooth codes:
*  m = monotonic
*  p = p-spline
*  i/d = increasing/decreasing
*  cv = concavity

unimodal splines? (uniReg package)
