setwd("~/Projects/UiS-Git/it-ift/implementation")
library(TMB)
compile("jd_ift_exact_spa.cpp")
dyn.load(dynlib("jd_ift_exact_spa"))

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
    data <- list(X=log(X), dt=1/12, process=2, scheme=1, jump=0, niter=3, ghiter=30)
    plot(X, type="l", main="Simulated GBM")
}

par_diff <- c(0.2,0.1)
par_jump <- c()
param <- list(par = c(par_diff,par_jump))
obj <- MakeADFun(data, param)
obj$fn(obj$par)
opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he, control=list(trace=1))
obj$fn(opt$par)
res <- sdreport(obj)
res


dyn.unload(dynlib("jd_ift_exact_spa"))
