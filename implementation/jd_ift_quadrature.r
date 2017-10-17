setwd("~/Projects/UiS-Git/it-ift/implementation")
library(TMB)
compile("jd_ift_quadrature.cpp")
dyn.load(dynlib("jd_ift_quadrature"))

# real data
real_data <- TRUE
if(real_data){
    library(Quandl)
    start_date <- "2013-01-01"; end_training <- "2015-01-01"; 
    test <- Quandl("BCB/UDJIAD1",trim_start=start_date)
    library(Quandl)
    DJIA<-Quandl("BCB/UDJIAD1",trim_start=start_date, trim_end=end_training)
    DJIA <- DJIA[rev(rownames(DJIA)),]
    plot(DJIA,type="l")
    log_price <- log(DJIA$Value)
    data <- list(X=log_price, dt=1/250, process=2, scheme=1, jump=0, qiter=100, quadrature=2)
}

# Simulated data
simulate = FALSE
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
    data <- list(X=log(X), dt=1/12, process=2, scheme=1, jump=0, qiter=100, quadrature=2)
    plot(X, type="l", main="Simulated GBM")
}

# GBM Estimation
gbm_estimation <- FALSE
if(gbm_estimation){
    par_diff <- c(kappa=0.1,sigma=0.2)
    par_jump <- c()
    param <- list(par = c(par_diff,par_jump))
    obj <- MakeADFun(data, param)
    opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he, control=list(trace=1))
    res <- sdreport(obj)
    res    
}


# MJD Estimation
mjd_estimation <- TRUE
if(mjd_estimation){
    data$process = 3
    data$jump = 1
    data$qiter = 180
    par_diff <- c(kappa=0.1,sigma=0.2)
    par_jump <- c(30,0,0.05)
    param <- list(par = c(par_diff, par_jump))
    obj <- MakeADFun(data, param)
    obj$fn(obj$par)
    opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he, control=list(trace=1))
    res <- sdreport(obj)
    res    
}


dyn.unload(dynlib("jd_ift_quadrature"))