# Minimal it-ift 
# R part of TMB program
# Author: Berent Å. S. Lunde
# 19.09.2017
setwd("~/Projects/UiS-Git/it-ift/implementation")

library(TMB)

# Compile and load into R
compile("generic_itift_tmp.cpp")
dyn.load(dynlib("generic_itift_tmp"))

# Simulate data
simulate = TRUE
GBM <- TRUE
OU <- FALSE
if(simulate){
    if(GBM){
        source("simulation/Simulation_GBM.R")
        set.seed(123)
        time = 10
        N=250*time
        mu = 0.1
        sigma = 0.05
        x0 = 1
        dt <- time / N
        seed = 123
        X <- log(GBM_process(time, N, mu, sigma, x0, seed))
        plot(X, type="l", main="Simulated GBM")
        param <<- list( par = c(mu, sigma) )
        process <- 2
    }
    if(OU){
        source("simulation/Simulation_OU.R")
        set.seed(123)
        time = 30
        N=250*time
        dt = time / N
        kappa=2
        alpha=1
        sigma=1
        x0=3
        seed=123
        X <- OU_process(time,N,kappa,alpha,sigma,x0,seed)
        plot(X,type="l", main="Simulated OU")
        param <- list(par=c(kappa,alpha,sigma))
        process=3
    }
}

# Construct TMB object
data <- list(Xt=X, dt=dt, process=process, scheme=1, jump=0, ghiter=30)
#param <- list( par = c(0.2, 0.2) )
obj <- MakeADFun(data, param)

# test
obj$fn(obj$par)
obj$gr(obj$par)
obj$he(obj$par)

opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he, control = list(trace=1))
rep <- sdreport(obj)
rep

# Comparison
nll_r <- function(par, X, dt){
    nll = 0
    for(i in 2:length(X)){
        nll = nll - dnorm(X[i], X[i-1]+(par[1]-.5*par[2]^2)*dt, par[2]*sqrt(dt),log=TRUE)
    }
    return(nll)
}
nll_r(c(mu,sigma), X, dt)
opt_r <- nlminb(c(mu,sigma),nll_r, X=X, dt=dt)
opt_r

# Optimise

# Unload from R
dyn.unload(dynlib("generic_itift_tmp"))