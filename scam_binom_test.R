## make sure we have the latest/bug-fixed versions
remotes::install_github("bbolker/reformulas")
remotes::install_github("glmmTMB/glmmTMB/glmmTMB")

library(glmmTMB)
library(scam)
set.seed(101)
dd <- data.frame(x=seq(-5, 5, length = 101))
dd <- within(dd, {
             p <- plogis(1 + x - x^3/4)
             b0 <- rbinom(nrow(dd), size = 1, prob = p)
             b1 <- rbinom(nrow(dd), size = 20, prob = p)
             })

par(las = 1, bty = "l") ## cosmetic
colvec <- c(2,4:6, 8)

## all approaches behave similarly for Bernoulli data
## scam only does GCV, glmmTMB only does ML/REML
m_bern_gam_gcv <- gam(b0 ~ s(x, bs = "tp"), family = binomial, data = dd,
                      method = "GCV.Cp")
m_bern_gam_reml <- gam(b0 ~ s(x, bs = "tp"), family = binomial, data = dd,
                       method = "REML")
m_bern_scam <- scam(b0 ~ s(x, bs = "tp"), family = binomial, data = dd)
m_bern_glmmTMB <- glmmTMB(b0 ~ s(x, bs = "tp"), family = binomial, data = dd,
                          REML = TRUE)
predmat_bern <- cbind(predict(m_bern_gam_gcv),
                 predict(m_bern_gam_reml),
                 predict(m_bern_scam),
                 predict(m_bern_glmmTMB))
matplot(dd$x, predmat_bern, lty = 1:4, col = colvec)
legend(lty = 1, lwd = 2,
       col = colvec[1:4],
       x = 1.5, y = 8,
       legend = c("gam/GCV", "gam/REML", "scam/GCV", "glmmTMB/REML"))

lines(dd$x, qlogis(dd$p), lwd = 2)
rug(dd$x[dd$b0==0], ticksize = 0.1, side = 1)
rug(dd$x[dd$b0==1], ticksize = 0.1, side = 3)
## results are heavily smoothed relative to true prob, but not surprising given
## v. low-resolution data

m_binom_gam_gcv <- gam(cbind(b1, 20-b1) ~ s(x, bs = "tp"), family = binomial, data = dd,
                      method = "GCV.Cp")
m_binom_gam_reml <- gam(cbind(b1, 20-b1) ~ s(x, bs = "tp"), family = binomial, data = dd,
                       method = "REML")
m_binom_scam_2col <- scam(cbind(b1, 20-b1) ~ s(x, bs = "tp"), family = binomial, data = dd)
## NOTE warnings() about non-integer # successes (due to weights issues)
m_binom_scam_wts <- scam(b1/20 ~ s(x, bs = "tp"), weights = rep(20, nrow(dd)),
                                                             family = binomial, data = dd)
m_binom_glmmTMB <- glmmTMB(cbind(b1, 20-b1) ~ s(x, bs = "tp"), family = binomial, data = dd,
                           REML = TRUE)

expand_bern <- function(dd, response = "b1", size = 20) {
    apply(dd, 1,
          function(x) {
              browser()
              k <- x[[response]]
                       data.frame(response = rep(c(0,1), times = c(size-k, k)),
                                  x[names(x) != "response"])
          })
}

m_binom_scam_expand <- scam(b1/20 ~ s(x, bs = "tp"), weights = rep(20, nrow(dd)),
                                                             family = binomial, data = dd)



predmat_binom <- cbind(predict(m_binom_gam_gcv),
                 predict(m_binom_gam_reml),
                 predict(m_binom_scam_2col),
                 predict(m_binom_scam_wts),
                 predict(m_binom_glmmTMB))



matplot(dd$x, predmat_binom, lty = 1:4, col = colvec, ylim = c(-25,25))
lines(dd$x, qlogis(dd$p), lwd = 2)
## squeeze obs probs in slightly so we can look at logit scale
points(dd$x, qlogis((dd$b1+0.25)/20.5))
legend("topright", lty = 1, lwd = 2,
       col = colvec,
       legend = c("gam/GCV", "gam/REML", "scam/GCV/2col",
                  "scam/GCV/2col/wts", "glmmTMB/REML"))



