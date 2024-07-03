expand_bern <- function(dd, response = "b1", size = 20) {
    dd[["..size"]] <- if (is.numeric(size)) size else dd[[size]]
    L <- apply(dd, 1,
               function(x) {
                   k <- x[[response]]
                   data.frame(response = rep(c(0,1), times = c(x[["..size"]]-k, k)),
                              rbind(x[!names(x) %in% c(response, "..size")]))
          })
    ret <- do.call(rbind, L)
    names(ret)[names(ret)=="response"] <- response
    return(ret)
}

holling2_rtmb <- function(parms) {
    getAll(tmbdat, parms)
    prob <- a/(1+a*h*Initial)
    logitprob <- 1/(1+exp(-prob))
    Killed %~% dbinom(prob = prob, size = Initial)
    REPORT(prob)
    ADREPORT(logitprob)
}

## parameters <- list(a = 1, h = 50)
## tmbdat <- dd
## obj <- MakeADFun(holling2_rtmb, parameters)
## obj$fn()
## res <- with(obj, nlminb(par, fn, gr), control = list(eval.max = 500))
## obj$report()$prob
## sdreport(obj)$value
## sdreport(obj)$sd
