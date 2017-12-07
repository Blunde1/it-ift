# R program that builds the likelihood of a jump diffusion based upon
# exact spa methodlolgy.
# Note: for real data estimation - build independent R script
# Part 1. Simulation
# Part 2. Likelihood estimation
# Part 3. Transition density plotting
# Part 4. Sampled likelihood estimates

setwd("C:/Users/Berent/Projects/it-ift/implementation")
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
        par_diff <- c(mu+0.5,sigma/10)
        par_jump <- numeric(0)
        param <- list(par=par_diff)
        data <- list(X=log(X), dt=dt, process=2, scheme=1, jump=0, niter=3, ghiter=32, line_search=0)
        plot(X, type="l", main="Simulated GBM")
    }
    if(mjd){
        source("simulation/Simulation_MJD.R")
        set.seed(123)
        time = 2
        N=250*time
        r = 0.55
        sigma = 0.2
        lambda = 30
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
        data <- list(X=log(X$x), dt=dt, process=3, scheme=1, jump=1, niter=50, ghiter=180, line_search=0)
    }
}

##### REAL DATA ####
library(Quandl)
start_date <- "1950-01-01"; end_training <- "2017-01-01"; 
test <- Quandl("BCB/UDJIAD1",trim_start=start_date)
library(Quandl)
DJIA<-Quandl("BCB/UDJIAD1",trim_start=start_date, trim_end=end_training)
#DJIA <- Quandl("BCB/UDJIAD1")
DJIA <- DJIA[rev(rownames(DJIA)),]
plot(DJIA,type="l")
log_price <- log(DJIA$Value)
data <- list(X=log_price, dt=1/250, process=3, scheme=1, jump=1, niter = 5, qhiter = 16, line_search=0)

#####
# Part 2. Likelihood estimation
#par_diff <- c(0.1,0.2)
#par_jump <- c(30,0,0.05)
#param <- list(par = c(par_diff,par_jump))
# data and param are set during simulation
# Works for real data!
obj <- MakeADFun(data, param)
obj$fn(obj$par)
obj$gr(obj$par)
obj$he(obj$par)
opt <- nlminb(opt$par, obj$fn, obj$gr, control=list(trace=1))
opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he, control=list(trace=1))
obj$fn(opt$par)
res <- sdreport(obj)
res

#####
# Part 3. Transition density plotting
td_settings <- TRUE
if(td_settings){
    r = 0.2
    log.sigma = log(0.3)
    log.lambda = log(100)
    mu=-0.01
    log.nu=log(0.05)
    x0 = 0
    dt <- 1/250
    par_diff <- c(r,log.sigma)
    par_jump <- c(log.lambda,mu,log.nu)
    param <- list(par=c(par_diff, par_jump))
    data <- list(X=numeric(2), dt=dt, process=3, scheme=1, jump=1, niter=30, ghiter=32, line_search = 0)
}
# Settings from master thesis
param$par <- c(0.4,log(0.3),log(10),-0.01,log(0.05)) # Actual setting
data$dt <- 1/250
# param$par <- c(0.55,log(0.2),log(18),-0.0063,log(0.4))

param <- list(par=c("r"=0.03,"sigma"=log(0.2),"lambda"=log(1),"mu"=-0.5,"nu"=log(0.1)))
data$dt <- 1/4
x_vals <- seq(-1.7,0.5,by=0.01)

x_vals <- seq(-0.12, 0.12, by=0.001)
y_vals <- numeric(length(x_vals))
for(i in 1:length(x_vals)){
    data$X <- c(x0,x_vals[i])
    obj <- MakeADFun(data, param)
    y_vals[i] <- exp(-obj$fn()) # obj returns nll
}
lines(x_vals,y_vals,col="green")
plot(x_vals, y_vals, type="l") # bimodal works

dyn.unload(dynlib("jd_ift_exact_spa"))
