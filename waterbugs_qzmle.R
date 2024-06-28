remotes::install_github("queezzz/qzmle")
x <- read.csv("Response_surfaces_data_Gamboa_data.csv")
x$block <- factor(x$block) #Make the block labels into factors
bugdat <-  subset(x,predtype=="belo", select=-predtype)  #Isolate data for Odonates

library(ggplot2); theme_set(theme_bw())
ggplot(bugdat, aes(initial,killed/initial,colour=block)) +
    stat_sum(alpha=0.5) +
    facet_wrap(~sizeclass) +
    scale_size(breaks=c(1,2),range=c(2,5))

library(qzmle)
library(bbmle)

form <- killed ~ dbinom(
            prob=a*exp(-size/b_size)/(1+a*exp(-size/b_size)*h*initial),
            size=initial)

system.time(m1 <- qzmle::mle(form,
                      data=bugdat,
                      start=list(a=1,b_size=1,h=1)
                      ))

system.time(m4 <- qzmle::mle(form,
                      data=bugdat,
                      start=list(a=1,b_size=1,h=1),
                      links=c(a="logit",h="log")
                      ))

system.time(m5 <- qzmle::mle(form,
                      data=bugdat,
                      start=as.list(coef(m2))
                      ))

system.time(m6 <- qzmle::mle(form,
                      data=bugdat,
                      start=list(a=1,b_size=1,h=1),
                      links=c(a="logit",h="log"),
                      method="TMB"
                      ))


## hmm, bbmle is faster ???
system.time(m2 <- bbmle::mle2(form,
                   data=bugdat,
                   control=list(maxit=5000),
                       start=list(a=1,b_size=1,h=1)
                   ))


## with Nelder-Mead
system.time(m3 <- bbmle::mle2(form,
                   size=initial),
                   data=bugdat,
                   method="Nelder-Mead",
                   control=list(maxit=5000),
                       start=list(a=1,b_size=1,h=1)
                   ))
