## testing glmmTMB etc. with simple mpd examples

library(scam)
library(glmmTMB)
set.seed(3)
n <- 100
x <- runif(n)*3-1
f <- exp(-1.3*x)
y <- rpois(n,exp(f))
dat <- data.frame(x=x,y=y)
## fit model ...
b0 <- scam(y~s(x,k=15,bs="mpd"),family=poisson(link="log"),
            data=dat)
## unconstrained model fit for comparison...
b1 <- scam(y~s(x,k=15,bs="ps"),family=poisson(link="log"),
              data=dat)

b3 <- glmmTMB(y~s(x,k=15,bs="mpd"),family=poisson(link="log"),
            data=dat)
## unconstrained model fit for comparison...
b4 <- glmmTMB(y~s(x,k=15,bs="ps"),family=poisson(link="log"),
              data=dat)

xvec <- seq(-1, 2, length = 101)
pf <- data.frame(x = xvec)

## fails
plot(xvec, predict(b3, newdata = pf))



## first smooth term of first RE of each model
ss3 <- b3$modelInfo$reTrms[[1]]$smooth_info[[1]]
ss4 <- b4$modelInfo$reTrms[[1]]$smooth_info[[1]]

devtools::load_all("~/R/pkgs/glmmTMB/glmmTMB")
debug(getXReTrms)
plot(xvec, predict(b4, newdata = pf))

fr <- model.frame(b4)
fr <- model.frame(~x, data = pf)
ff(ss4)

xr3 <- ff(ss3)
beta <- fixef(b3)$cond
b <- getME(b3, "b")
pred <- beta[1] + xr3$re$Xf*beta[2] + xr3$re$rand$Xr %*% b
plot(xvec, pred)
lines(xvec, predict(b4, newdata = pf))


## function from getXReTrms
ff <- function(s) {
    if (is.null(s)) return(NULL)
    X <- PredictMat(s$sm, fr)   ## get prediction matrix for new data
    if (all(abs(X[,1] - 1.0) <1e-12)) X <- X[,-1] ## drop intercept
    ## transform to r.e. parameterization
    if (!is.null(s$re$trans.U)) X <- X%*%s$re$trans.U
    X <- t(t(X)*s$re$trans.D)
    ## re-order columns according to random effect re-ordering...
    X[,s$re$rind] <- X[,s$re$pen.ind!=0] 
    ## re-order penalization index in same way  
    pen.ind <- s$re$pen.ind; s$pen.ind[s$re$rind] <- pen.ind[pen.ind>0]
    ## start return object...
    s_new <- list(re = list(rand=list(), Xf=X[,which(s$re$pen.ind==0),drop=FALSE]))
    for (i in 1:length(s$re$rand)) { ## loop over random effect matrices
        s_new$re$rand[[i]] <- X[, which(pen.ind==i), drop=FALSE]
        attr(s_new$re$rand[[i]], "s.label") <- attr(s$re$rand[[i]], "s.label")
    }
    names(s_new$re$rand) <- names(s$re$rand)
    return(s_new)
}


ff(ss3)

## how is constructed smooth different
