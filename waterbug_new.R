library(tidyverse); theme_set(theme_bw())
library(plotly)
library(bbmle)
library(emdbook)
library(mgcv)
library(scam)

source("funs.R")
## rayshader is disappointing/frustrating
## library(rayshader)

datfn <- "McCoy_response_surfaces_Gamboa.csv"

x <- (read.csv(datfn)
    |> transform(block = factor(block))
    |> subset(cohort == "single" & predtype == "belo",
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
    geom_point(aes(colour = killed/initial), size = 4)
    ## scale_colour_viridis_c()
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

p2 <- (plot_ly(x= ~initial, y = ~size, z = ~prop)
    |> add_markers(data = x, marker = marker)
    |> add_paths(data = seg_data(x, prop))
    |> layout(scene = list(yaxis = list(rangemode = "tozero"),
              xaxis = list(rangemode = "tozero")))
)
print(p2)

## https://www.datanovia.com/en/blog/how-to-create-a-ggplot-like-3d-scatter-plot-using-plotly/
## https://community.plotly.com/t/droplines-from-points-in-3d-scatterplot/4113/10

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

p2 |> add_trace(type =  "mesh3d", data = long_fmt(cc), opacity = 0.4)

fit2 <- gam(cbind(killed, initial-killed) ~ te(size, initial),
            data = x, family = binomial)

gam_pred <- expand.grid(size = 0:60, initial = 0:100)
gam_pred$prop <- predict(fit2, newdata = gam_pred, type = "response")

p2 |> add_trace(type =  "mesh3d", data = gam_pred, opacity = 0.4)

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

