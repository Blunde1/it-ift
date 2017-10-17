# R program that builds the likelihood of a jump diffusion based upon
# exact spa methodlolgy.
# Note: for real data estimation - build independent R script
# Part 1. Simulation
# Part 2. Likelihood estimation
# Part 3. Transition density plotting
# Part 4. Sampled likelihood estimates


setwd("~/Projects/UiS-Git/it-ift/implementation")
library(TMB)
compile("jd_ift_exact_spa.cpp")
dyn.load(dynlib("jd_ift_exact_spa"))

##### 
# Part 1 - SIMULATION
# Simulated data
simulate = TRUE
gbm <- FALSE
mjd <- TRUE
if(simulate){
    if(gbm){
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
        par_diff <- c(mu,sigma)
        par_jump <- numeric(0)
        data <- list(X=log(X), dt=1/12, process=2, scheme=1, jump=0, niter=3, ghiter=30)
        plot(X, type="l", main="Simulated GBM")
    }
    if(mjd){
        source("simulation/Simulation_MJD.R")
        set.seed(123)
        time = 2
        N=250*time
        r = 0.55
        sigma = 0.2
        lambda = 18
        mu=-0.01
        nu=0.04
        x0 = 1
        dt <- time / N
        seed = 123
        X <- mjd_process(N,time,r,sigma,lambda,mu,nu,x0,seed)
        par_true <- c(r,sigma,lambda,mu,nu)
        par_diff <- c(r,sigma)
        par_jump <- c(lambda,mu,nu)
        param <- list(par=c(par_diff, par_jump))
        data <- list(X=log(X$x), dt=dt, process=3, scheme=1, jump=1, niter=50, ghiter=32,alpha=1)
    }
}

#####
# Part 2. Likelihood estimation
#par_diff <- c(0.1,0.2)
#par_jump <- c(30,0,0.05)
#param <- list(par = c(par_diff,par_jump))
# data and param are set during simulation
obj <- MakeADFun(data, param)
obj$fn(obj$par)
opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he, control=list(trace=1))
obj$fn(opt$par)
res <- sdreport(obj)
res

#####
# Part 3. Transition density plotting
td_settings <- TRUE
if(td_settings){
    r = 0.4
    sigma = 0.3
    lambda = 10
    mu=-0.01
    nu=0.05
    x0 = 0
    dt <- 1/250
    par_diff <- c(r,sigma)
    par_jump <- c(lambda,mu,nu)
    param <- list(par=c(par_diff, par_jump))
    data <- list(X=numeric(2), dt=dt, process=3, scheme=1, jump=1, niter=50, ghiter=180,alpha=1)
}
x_vals <- seq(-0.09, 0.09, by=0.001)
y_vals <- numeric(length(x_vals))
for(i in 1:length(x_vals)){
    data$X <- c(x0,x_vals[i])
    obj <- MakeADFun(data, param)
    y_vals[i] <- exp(-obj$fn()) # obj returns nll
}
lines(x_vals,y_vals,col="green")
plot(x_vals, y_vals, type="l") # bimodal works

dyn.unload(dynlib("jd_ift_exact_spa"))
