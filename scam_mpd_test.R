remotes::install_github("bbolker/reformulas")
remotes::install_github("glmmTMB/glmmTMB/glmmTMB")

library(glmmTMB)
library(mgcv)
library(scam)
library(Matrix)
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

sm <- smoothCon(s(x, bs = "mpd"), data = dd, absorb.cons = TRUE)
names(sm[[1]])

image(Matrix(sm[[1]]$S[[1]]))
image(Matrix(sm[[1]]$P[[1]]))
