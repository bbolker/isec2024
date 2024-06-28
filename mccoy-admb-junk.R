do_admb("mccoypredb",inp[-3],
        list(a=0.1,h=1))

do_admb("simple",list(nobs=nrow(x),Y=x$killed,x=x$initial),
        list(a=0.1,b=1.05))

do_admb("mccoypredc",list(nobs=nrow(x),Y=x$killed,x=x$initial),
        list(h=0.1,a=1.05))
do_admb("mccoypredd",list(nobs=nrow(x),Y=x$killed,x=x$initial),
        list(h=0.1,a=1.05))

coef(lm(killed~initial,data=x))
do_admb("mccoypredb",inp[-3],
        list(h=0.25,a=1.05))

coef(lm(killed~initial+I(initial^2)-1,data=x))
coef(lm(killed~initial+I(initial^2)-1,data=ododat))
do_admb("mccoyprede",inp[-3],
        list(h=-0.0011,a=0.317))

do_admb("mccoypredf",inp[-3],
        list(h=-0.002,a=0.32))

do_admb("mccoypredf",inp[-3],
        list(h=-0.0025,a=0.47))

file.copy("mccoypredf.dat","mccoypredf.dat1")
do_admb("mccoypredf",list(nobs=nrow(x),killed=x$killed,initial=x$initial),
        list(h=-0.002,a=0.32))

## works/fails depending on whether we put in the whole data set or not

do_admb("mccoypredg",list(nobs=nrow(x),killed=x$killed,initial=x$initial,
                          init2=x$initial^2),
        list(h=-0.002,a=0.32))


set.seed(1)
initial = rep(10,10)
killed = rbinom(10,prob=0.4,size=initial)
  
do_admb("mccoypreda",list(nobs=10,killed=killed,initial=initial),
        list(p=0.3))
do_admb("mccoypreda",list(nobs=nrow(x),killed=x$killed,initial=x$initial),
        list(p=0.3))
sum(x$killed)/sum(x$initial)

-sum(dbinom(killed,size=initial,prob=0.3,log=TRUE))
-sum(dbinom(killed,size=initial,prob=0.4,log=TRUE))
-sum(killed*log(0.3)+(initial-killed)*log(1-0.3))
-sum(killed*log(0.4)+(initial-killed)*log(1-0.4))
## seems well-behaved
prob = with(c(pars,as.list(ododat)),1/(1/(c*exp(-size/d))+h*initial))
summary(prob)
binloglik = with(ododat,lchoose(initial,killed)+killed*log(prob)+
  (initial-killed)*log(1-prob))
summary(binloglik)
sum(binloglik)
par_read("mccoypred0")




for (i=Z.rowmin(); i<=Z.rowmax(); i++) {
     cvec(i) = c;
     for (j=Z.colmin(); j<=Z.colmax(); j++) {
      cvec(i) += sigma_c*Z(i,j)*u(j);
     }
  }
  //  cvec = sigma_c*cvec + c;

  cout << "test calc \n";
  // u + u;
  cout << "test calc 1A \n";
  // cout << u + 1;
  cout << "test 2 \n";
  // cout << 
  Z*u;
  cout << "test 3 \n";
  cvec = sigma_c*Z*u;
  cout << cvec.indexmin();
  cout << "\ntest 4 \n"; 
  cout << "cvec calc \n";
  
      // cout << "v(i): ";
      // cout << v(i);
      // cout << "\nx: ";
      // cout << x;
      // cout << "\nres: ";
      // cout << tmp(i);
      // cout << "\n";




\begin{verbatim}
export ADMBSRCDIR=~/lib/admb-project-read-only/
export PROJDIR=~/students/mccoy
export ADMB_HOME=/usr/local/src/admb
adcomp -r f1b2vc1
## g++ -s  -L/usr/local/src/admb/lib mccoypred5.o \
##    -ldf1b2s -ladmod -ladt -lado -ldf1b2s -ladmod -ladt -lado -o mccoypred5
cd $ADMBSRCDIR/df1b2-separable
export CFLAGS="-O3 -Wno-deprecated -D__GNUDOS__  -Dlinux -DOPT_LIB -DUSE_LAPLACE -fpermissive -I. -I$ADMB_HOME/include"
ar -rs libdf1b2o.a *.o

g++ -c $CFLAGS  *.cpp  ## takes a while
cp libdf1b2o.a ~/students/mccoy
g++ -s  -L/usr/local/src/admb/lib mccoypred5.o -ladmod -ladt -lado  -ldf1b2s -o mccoypred5w

/usr/local/src/admb/lib/libdf1b2s.a(f1b2vc1.obj): In function `df1b2vector::df1b2vector()':
f1b2vc1.cpp:(.text+0x0): multiple definition of `df1b2vector::df1b2vector()'
mccoypred5.o:mccoypred5.cpp:(.text+0x20): first defined here

cp /usr/local/src/admb/lib/libdf1b2s.a .
ar -t libdf1b2s.a

g++ -s  -L. -I$ADMBSRCDIR/nh99 -I$ADMBSRCDIR/linad99 -I$ADMBSRCDIR/tools99 -L/usr/local/src/admb/lib mccoypred5.o -ladmod -ladt -lado  -ldf1b2s_hack -o mccoypred5
\end{verbatim}

\item{going a little nuts making up variants of things trying
    to find out what's wrong with my code.  Probably something
    really boneheaded;
     \begin{itemize}
     \item \code{mccoypred0}: ``simple'' version of model,
       no random effects.  Have simplified as far as an
       exponential dependence of attack rate on size.
       FAILS.
     \item \code{mccoypreda}: binomial likelihood with
       a constant probability of mortality (size, density-ind.).
       WORKS.
     \item \code{simple}: regression code from ADMB examples.
       WORKS with initial as a linear function of killed.
     \item \code{mccoypredb}: attempted regression model of my
       own. FAILS.
     \item \code{mccoypredc}: same as mccoypredb, with
       variable names changed from killed, initial, pred_killed
       to x, Y, pred_Y. WORKS. (!!!!????)
     \item \code{mccoypredd}: change x back to initial (WORKS);
       change pred_Y back to pred_killed (WORKS);             
     \item \ldots !!!! I screwed up the order of the 
       parameters again \ldots \code{mccoypredb} works with
       the parameters in the right order.
     \item \code{mccoyprede}: attempt to define a density-dependent
       probability and multiply it elementwise by pred\_killed. FAILS.
     \item \code{mccoypredf}: attempt to do the same thing as
       a quadratic regression, using \code{square} or \code{elem\_prod(initial,initial)}. FAILS.
     \item \code{mccoypredg}: ditto, with density-squared computed in
       R. WORKS.
       \code{mccoypredb}
     \end{itemize}}

\section{More pictures}
Pictures for talk: (predict does not work when parameters
specified!)

<<pred3>>=
sizemeans = tapply(x$size,list(x$sizeclass),mean)
predframe = expand.grid(predtype=factor(c("belo","odo")),
  sizeclass=1:5,initial=seq(2,90,by=2),block=factor(1))
predframe$size = sizemeans[predframe$sizeclass]
pred1 = predict(g3_exp,newdata=predframe,type="response")
predframe = data.frame(predframe,prob=pred1,killed=pred1)
@ 

Dashed lines are
smoothed nonparametric curves (conf intervals omitted
for clarity), solid lines are glm fit.

<<plot_pred3,fig=TRUE>>=
belodat =  subset(x,predtype=="belo")
## x2 = ododat
x2 = belodat
x2$sizeclass = factor(x2$sizeclass)
g3_expA = mle2(killed~dbinom(prob=1/(1/(c*exp(-size/d))+
                                 h*initial),
  size=initial),
  start=sv,
  data=x2)

p1 = ggplot(x2,aes(x=initial,y=killed,group=sizeclass,
  fill=sizeclass,
  colour=sizeclass))+ss
print(p1+geom_point(aes(pch=block))+
                    facet_wrap(~predtype)+
                    stat_sum(aes(size=..n..))+
                    geom_smooth())

sizemeans = tapply(x2$size,list(x2$sizeclass),mean)
predframe = expand.grid(sizeclass=levels(x2$sizeclass),
  initial=seq(0,100,by=2),block=factor(1))
predframe$size = sizemeans[as.numeric(predframe$sizeclass)]
pred1 = predict(g3_expA,newdata=predframe)
predframe = data.frame(predframe,prob=pred1,killed=pred1)

p2 = ggplot(x2,aes(x=size,y=killed,group=factor(initial),
  fill=factor(initial),
  colour=factor(initial)))+ss
@ 

<<fig=TRUE>>=
print(p1+geom_point()+stat_sum(aes(size=..n..))+
      geom_line(data=predframe))
@ 

<<fig=TRUE>>=
predframe = expand.grid(size=5:55,
  initial=unique(x2$initial),block=factor(1))
pred1 = predict(g3_expA,newdata=predframe)
predframe = data.frame(predframe,prob=pred1,killed=pred1)

print(p2+geom_point()+
      geom_line(data=predframe))
@ 
