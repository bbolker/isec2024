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
}
