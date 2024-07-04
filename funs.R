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

fit_mpd_fun <- function(data, response = "y",
                        xvar = "x",
                        form = s(x, bs = "mpd"), 
                        size = numeric(0),
                        parms = NULL,
                        smoothdata = data,
                        predict = FALSE,
                        family = "gaussian", random = "b1", silent = TRUE,
                        ...) {
    form$term <- xvar
    ## if predicting, make sure to set up smooth with old data
    sm1 <- smoothCon(form, data = smoothdata, absorb.cons = TRUE)[[1]]
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
    if (predict) {
        ## shouldn't need to map() b if random = NULL ?
        obj$fn(unlist(parms))
        sdr <- sdreport(obj)
        return(with(sdr,
                    data.frame(nm = names(value), value, sd)))
    }
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

mk_holling2_rtmb_fun <- function(data) {
    f <- function(parms) {
        getAll(data, parms)
        a <- exp(loga)
        prob <- a/(1+a*exp(logh)*Initial)
        logitprob <- qlogis(prob) ## for CIs on logit scale
        nll <- -1*sum(dbinom(Killed, prob = prob, size = Initial, log = TRUE))
        REPORT(prob)
        ADREPORT(logitprob)
        return(nll)
    }
    return(f)
}

fit_RTMB_holling2 <- function(data, parms = list(loga=log(0.5), logh = log(0.01))) {
    obj <- MakeADFun(mk_holling2_rtmb_fun(data), parameters = parms,
                     silent = TRUE)
    fit <- with(obj, nlminb(par, fn, gr))
    sdr <- sdreport(obj)
    return(list(obj = obj, fit = fit, pred = data.frame(Initial = data$Initial,
                                                        prob = obj$report()$prob,
                                                        logitprob = sdr$value,
                                                        logitprob_sd = sdr$sd)))
}

pred_RTMB_holling2 <- function(data, newdata, parms, response = "Killed") {
    newdata[[response]] <- NA_integer_
    n_new <- nrow(newdata)
    data <- rbind(newdata, data)
    ## don't think I have to worry about mapping b parameters here ...
    ## no random coef here
    newobj <- MakeADFun(mk_holling2_rtmb_fun(data), parameters = as.list(parms), silent = TRUE)
    newobj$fn(parms)
    sdr <- sdreport(newobj)    
    pframe <- with(sdr,
                   data.frame(Initial = data$Initial,
                              prob = newobj$report()$prob,
                              logitprob = value,
                              logitprob_sd = sd,
                              lwr = plogis(value-2*sd),
                              upr = plogis(value+2*sd))
                   )
    pframe <- pframe[seq(n_new), ]
    return(pframe)
}
    

