library(tidyverse); theme_set(theme_bw())
library(plotly)
library(bbmle)
library(emdbook)
library(mgcv)
library(scam)
if (packageVersion("scam") < "1.2.17.9000") stop("need up-to-date/hacked version, see BMB github")
library(RTMB)

source("funs.R")
## rayshader is disappointing/frustrating
## library(rayshader)

## other 3d options? rgl, scatterplot3d, ... ?

datfn <- "McCoy_response_surfaces_Gamboa.csv"

## belostomatids only
x <- (read.csv(datfn)
    |> transform(block = factor(block))
    |> subset(cohort == "single" & predtype == "odo",
              select = -c(cohort, predtype))
    |> droplevels()
    |> transform(csize = size - mean(size), prop = killed/initial)
)

## ar = f(c(block), size)
## h  = f(h*size)
## prob = f(ar, h)

gg1 <- ggplot(x, aes(x = initial, y = size)) +
    geom_point(aes(colour = killed)) +
    scale_colour_viridis_c()

gg2 <- ggplot(x, aes(x = initial, y = size)) +
    ## stat_sum(aes(colour= killed/initial)) +
    geom_point(aes(colour = killed/initial), size = 4) +
    scale_colour_viridis_c()
print(gg2)

marker <- list(color = ~prop,
               colorscale = c('#FFE1A1', '#683531'), 
               showscale = TRUE)

seg_data <- function(x, zvar) {
    ## why do we need enquo here??
    zvar <- enquo(zvar)
    xx <- (x
        |> mutate(.id = seq(nrow(x)))
        |> reframe(!!zvar := c(!!zvar, 0), across(-!!zvar), .by = .id)
        |> plotly::group2NA(".id")
    )
    return(xx)
}

## https://stackoverflow.com/questions/72281954/keep-other-columns-when-doing-group-by-summarise-with-dplyr
## https://stackoverflow.com/questions/50012328/r-plotly-showlegend-false-does-not-work
wb_plotly_0 <- (plot_ly(x= ~initial, y = ~size, z = ~prop)
    |> add_markers(data = x, marker = marker, showlegend = FALSE)
    |> add_paths(data = seg_data(x, prop), , showlegend = FALSE)
    |> hide_colorbar()
    |> layout(scene = list(yaxis = list(rangemode = "tozero"),
                           xaxis = list(rangemode = "tozero"),
                           camera = list(eye = list(x = 2.5, y = 2, z = 1)),
                           showlegend=FALSE))
)
if (interactive()) print(wb_plotly_0)

img <- function(obj, width = 1000, height = 1000, trim = TRUE) {
    nm <- paste0(deparse(substitute(obj)), ".png")
    save_image(obj, nm, width = width, height = height)
    if (trim) {
        system(sprintf("convert %s -trim tmp.png; mv tmp.png pix/%s", nm, nm))
    }
}
img(wb_plotly_0)

## https://www.datanovia.com/en/blog/how-to-create-a-ggplot-like-3d-scatter-plot-using-plotly/
## https://community.plotly.com/t/droplines-from-points-in-3d-scatterplot/4113/10

## are these for belo or odo??
L <- load("waterbug_fits_2.RData")
L <- load("waterbug_fits_2z.RData")
aictab <- tibble(resframe, aic = sapply(res, AIC)) |> mutate(across(aic, ~ . - min(.))) |> arrange(aic)
    

## power-ricker or  ricker attack rate, proportional handling time
fit1 <- mle2(killed ~ dbinom(prob = 1/(1/(c*size/d*exp(1-size/d)) + h*size*initial),
                             size = initial),
             start = list(c=1,d=20,h=20),
             control = list(parscale= c(c=1,d=20,h=20)),
             data = x,
             method = "Nelder-Mead")

long_fmt <- function(cc, nm = c("size", "initial", "prop")) {
    tibble(x = rep(cc$x, ncol(cc$z)),
           y = rep(cc$y, each = nrow(cc$z)),
           z = c(cc$z)) |>
        setNames(nm)
}
    
cc <- curve3d(1/(1/(c*size/d*exp(1-size/d)) + h*size*initial),
              xlim = c(0, 60),
              ylim = c(0, 100),
##            1/(1/(c*size/d*exp(1-size/d)) + h*size*initial)
        data = as.list(coef(fit1)),
        varnames = c("size", "initial"),
        sys3d = "image")


## https://stackoverflow.com/questions/34178381/how-to-specify-camera-perspective-of-3d-plotly-chart-in-r
wb_plotly_param <- (wb_plotly_0
    |> add_trace(type =  "mesh3d", data = long_fmt(cc), opacity = 0.4)
)
img(wb_plotly_param)

wb_gam_te <- gam(cbind(killed, initial-killed) ~ te(size, initial),
            data = x, family = binomial)

wb_gam_pred <- expand.grid(size = 0:60, initial = 0:100)
wb_gam_pred$prop <- predict(wb_gam_te, newdata = gam_pred, type = "response")

wb_plotly_gam <- (wb_plotly_0
    |> add_trace(type =  "mesh3d", data = wb_gam_pred, opacity = 0.4)
)
img(wb_plotly_gam)

xx <- x |> select(killed, size, initial) |> expand_bern(response = "killed", size = x$initial)
## tesmd2, tesmd1 = smooth monotone decreasing in var1/2, no constraint otherwise
## fit3 <- scam(killed ~ s(initial, size, bs = "tesmd2"), data = xx, family = binomial)
## fit3 <- scam(killed ~ s(initial, size, bs = "tesmd1"), data = xx, family = binomial)
## fails for belo .... (inner loop 3, can't correct step size)

## https://plotly.com/r/static-image-export/

## decreasing wrt var 1, convex wrt var 2
wb_scam_tedecv <- scam(killed ~ s(initial, size, bs = "tedecv"), data = xx, family = binomial)

scam_pred <- gam_pred
scam_pred$prop <- predict(fit3, newdata = scam_pred, type = "response")
wb_plotly_scam <- wb_plotly_0 |> add_trace(type =  "mesh3d", data = scam_pred, opacity = 0.4)
img(wb_plotly_scam)

ss <- s(initial, size, bs = "tedecv")
x2 <- x |> rename(Size = "size") ## hack, 'size' is confounded
fit_RTMB <- fit_mpd_fun(data = x2[c("killed", "Size", "initial")],
            response = "killed",
            xvar = c("Size", "initial"),
            form = s(size, initial, bs = "tesmd2"),
            size = x2$initial,
            family = "binomial",
            inner.control = list(smartsearch=FALSE, maxit =1),
            opt = "BFGS",
            start = list(b0 = -2, log_smD = 2),
            ## works with random = NULL; start from better values?
            random = NULL)

xp <- expand.grid(Size = 8:60, initial = 6:100)
k <- (smoothCon(s(Size, initial, bs="tesmd2"),
                data = x2, absorb.cons = TRUE)[[1]]$knots
    |> as.data.frame()
    |> setNames(c("Size", "initial"))
)

if (FALSE) {
    preds0 <- fit_mpd_fun(data = xp, response = "Killed",
                          size = x$initial, xvar = c("Size", "initial"),
                          family = "binomial", random = NULL,
                          knots = k, predict = TRUE,
                          parms = with(fit_RTMB$obj$env, parList(last.par.best)))
}

with(x, plot(fit_RTMB$mu, killed/initial))
abline(a=0, b=1)  ## could be worse?
## predictions?

## marginaleffects::predictions(fit3, newdata = data.frame(size=20, initial=20))
## devtools::load_all("~/R/pkgs/scam")
## predict(fit3, newdata = data.frame(size=20, initial=20))
## debug(predict.scam)

## X[, object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag %*% object$smooth[[k]]$Zc

## fp <- object$smooth[[k]]$first.para  ## 2
## lp <- object$smooth[[k]]$last.para ## 49
## X1 <- X[,fp:lp, drop = FALSE]
## object$smooth[[k]]$cmX

## X[, object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- sweep(X[, object$smooth[[k]]$first.para:object$smooth[[k]]$last.para], 2, object$smooth[[k]]$cmX)

cc_scam <- curve3d(predict(fit3, newdata = data.frame(size, initial), type = "response"),
              xlim = c(0, 60),
              ylim = c(0, 100),
##            1/(1/(c*size/d*exp(1-size/d)) + h*size*initial)
        varnames = c("size", "initial"),
        sys3d = "image")

if (interactive()) p2 |> add_trace(type =  "mesh3d", data = long_fmt(cc_scam), opacity = 0.4)

## unimodal??

## 
## for (i in 1:N) {
##     ## ar[i] <- cvec[block[i]]*pow(size[i]/d,gamma)*exp(1-size[i]/d)
##     ar[i] <- cvec[block[i]]*size[i]/d*exp(1-size[i]/d)
##     ## hfun[i] <- h*exp(hS*size[i])
##     hfun[i] <- h*csize[i]  ## prop.
##     prob[i] <- 1/(1/ar[i]+hfun[i]*initial[i])
##     killed[i] ~ dbin(prob[i],initial[i])
##   }
##   for (i in 1:nblock) {
##     cvec0[i] ~ dnorm(0,tau.c)
##     cvec[i] <- c*exp(cvec0[i])
##     killed.rep[i]~ dbin(prob[i],initial[i])
##   }
## ## priors
##   tau.c <- pow(sd.c,-2)
##   sd.c ~ dunif(0,1)
##   d ~ dlnorm(0,0.01)
##   gamma ~ dlnorm(0,0.01)
##   h ~ dunif(0,10) ## changed to positive uniform
##   ## hS ~ dunif(-1,1) 
##   c ~ dlnorm(0,0.01)

