library(tidyverse); theme_set(theme_bw())
library(rayshader)

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
plot_gg(gg1)

plot_gg(gg2)

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
seg_data(x, killed)

p <- plot_ly(xx, x= ~initial, y = ~size, z = ~killed, marker = marker) |>
    add_markers() 


p2 <- (plot_ly(x= ~initial, y = ~size, z = ~prop)
    |> add_markers(data = x, marker = marker)
    |> add_paths(data = seg_data(x, prop))
    |> layout(scene = list(yaxis = list(rangemode = "tozero"),
              xaxis = list(rangemode = "tozero")))
)
print(p2)

plot_ly(data = x, x = ~size, y = ~initial) |> add_markers() |>
    layout(yaxis = list(rangemode = "tozero"),
              xaxis = list(rangemode = "tozero"))
    
## segments?

## https://www.datanovia.com/en/blog/how-to-create-a-ggplot-like-3d-scatter-plot-using-plotly/
## https://community.plotly.com/t/droplines-from-points-in-3d-scatterplot/4113/10

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

