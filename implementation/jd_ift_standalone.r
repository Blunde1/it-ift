setwd("~/Projects/UiS-Git/it-ift/implementation")
library(TMB)
compile("jd_ift_standalone.cpp")
dyn.load(dynlib("jd_ift_standalone"))

# real data
real_data <- TRUE
if(real_data){
    start_date <- "2013-01-01"; end_training <- "2015-01-01"; 
    test <- Quandl("BCB/UDJIAD1",trim_start=start_date)
    library(Quandl)
    DJIA<-Quandl("BCB/UDJIAD1",trim_start=start_date, trim_end=end_training)
    DJIA <- DJIA[rev(rownames(DJIA)),]
    plot(DJIA,type="l")
    log_price <- log(DJIA$Value)
}

# Simulated data
simulate = TRUE
if(simulate){
    source("simulation/Simulation_GBM.R")
    set.seed(123)
    time = 100
    N=12*time
    mu = 0.1
    sigma = 0.2
    x0 = 1
    dt <- time / N
    seed = 123
    X <- GBM_process(time, N, mu, sigma, x0, seed)
    par_true <- c(mu,sigma)
    process=2; scheme=1; jump=0;
    plot(X, type="l", main="Simulated GBM")
}

param <- c(kappa=0.1,sigma=0.3)
parameters <- list(par=param)
data <- list(x=log(X), dt=1/12, process=1, scheme=1, jump=0, ghiter=100, quadrature=2)
obj <- MakeADFun(data, parameters)
opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he, control=list(trace=1))
res <- sdreport(obj)
res

# Comparison
nll_r <- function(par, X, dt){
    nll = 0
    for(i in 2:length(X)){
        nll = nll - dnorm(X[i], X[i-1]+(par[1]-.5*par[2]^2)*dt, par[2]*sqrt(dt),log=TRUE)
    }
    return(nll)
}
nll_r(c(mu,sigma), log(X), dt)
opt_r <- nlminb(c(mu,sigma),nll_r, X=log(X), dt=dt)
opt_r


#### Test with N(1,2) density - must comment out line 101 and comment in lines 102 and 103
f <- x <- seq(-5,8,length.out = 100)
for(i in 1:length(x)){
    data$x = c(x0,x[i])
    obj <- MakeADFun(data,parameters)
    f[i] <- obj$fn(obj$par)
}
plot(x,dnorm(x,1,2),type="l")
lines(x,exp(-f),type="l",col="red")


dyn.unload(dynlib("jd_ift_standalone"))