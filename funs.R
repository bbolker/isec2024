expand_bern <- function(dd, response = "y1", size = 20) {
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

## basic
mk_mpd_fun <- function(data, parms, random = "b1", silent = TRUE,
                       family = "gaussian", ...) {
    ## can't use %~% format if we want to add a penalty
    f <- function(parms) {
        getAll(data, parms)
        b_pos <- b1
        b_pos[p.ident] <- exp(b1)
        eta <- b0 + X %*% b_pos
        nll <- switch(family,
                      gaussian = { mu <- eta; -1*sum(dnorm(y, mu, exp(log_rSD), log = TRUE))},
                      binomial = { mu <- plogis(eta);
                          -1*sum(dbinom_robust(y, logit_p = eta, size = size, log = TRUE)) },
                      stop("unimplemented family"))
        ## translate from
        ## lambda = 1/sigma_sm^2
        ## MVgauss NLL = (1/2) (n*log(2pi) + log(det(Sigma)) + bT Sigma^{-1} b)
        ## Sigma = sd^2 S^{-1}
        ## log(det(Sigma)) = 2*logsd - log(det(S))
        ## MVNLL = C + logsd + 1/sd^2 (bT S b)
        pen <- (exp(-2*log_smSD) * (t(b_pos) %*% S %*% b_pos) + 2*log_smSD)/2
        REPORT(mu)
        ADREPORT(mu)
        REPORT(eta)
        ADREPORT(eta)
        nll + pen
    }
    MakeADFun(f, parms, random=random, silent = silent, ...)
}

fit_mpd_fun <- function(data, response = "y", form = s(x, bs = "mpd"), size = numeric(0), parms = NULL,
                        family = "gaussian", random = "b1", silent = TRUE, ...) {
    sm1 <- smoothCon(form, data = data, absorb.cons = TRUE)[[1]]
    parms <- parms %||% list(
                            b0 = 0,
                            b1 = rep(0, ncol(sm1$X)),
                            log_smSD = 0)
    if (family == "gaussian") parms <- c(parms, list(log_rSD = 0))
    data$y <- data[[response]]
    tmbdat <- c(as.list(data),
                list(size = size, family = family),
                list(p.ident = sm1$"p.ident", S = sm1$S[[1]], X = sm1$X))
    obj <- mk_mpd_fun(data = tmbdat, parms = parms,
                      random = random, 
                      silent = silent)
    ## p0 <- parms[names(parms) != "b1"]
    ## optim(par = p0, fn = obj$fn, control = list(maxit = 2000))
    res <- with(obj, nlminb(par, fn, gr))
    ret <- list(fit = res, obj = obj, mu = obj$report()$mu, eta = obj$report()$eta)
    class(ret) <- c("myRTMB", "list")
    return(ret)
}

predict.myRTMB <- function(x, type = "link", ...)  {
    switch(type,
           link = x$eta,
           response = x$mu)
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
