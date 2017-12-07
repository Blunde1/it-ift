# R program that builds the likelihood of a jump diffusion based upon
# exact mixed spa methodlolgy.
# Note: for real data estimation - build independent R script
# Part 1. Simulation
# Part 2. Likelihood estimation
# Part 3. Transition density plotting
# Part 4. Sampled likelihood estimates


setwd("C:/Users/Berent/Projects/it-ift/implementation")
setwd("~/Projects/UiS-Git/it-ift/implementation")
library(TMB)
compile("jd_exact_mspa.cpp")
dyn.load(dynlib("jd_exact_mspa"))


##### TRANSITION DENSITY ####
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
data <- list(X=numeric(2), dt=dt, process=3, scheme=1, jump=11, 
             niter_no_jump=5, niter_jump=64, ghiter_no_jump=16, ghiter_jump=32, 
             line_search = 0)

# Settings from master thesis
param$par <- c(0.4,log(0.3),log(10),-0.01,log(0.05)) # Actual setting
data$dt <- 1/250
# param$par <- c(0.55,log(0.2),log(18),-0.0063,log(0.4))

param <- list(par=c("r"=0.03,"sigma"=log(0.2),"lambda"=log(1),"mu"=-0.5,"nu"=log(0.1)))
data$dt <- 1/4
x_vals <- seq(-1.7,0.5,by=0.01)

x_vals <- seq(-0.10, 0.10, by=0.001)
y_vals <- numeric(length(x_vals))
for(i in 1:length(x_vals)){
    data$X <- c(x0,x_vals[i])
    obj <- MakeADFun(data, param)
    y_vals[i] <- exp(-obj$fn()) # obj returns nll
}
lines(x_vals,y_vals,col="green")
plot(x_vals, y_vals, type="l") # bimodal works

##### LIKELIHOOD ####

# Load data
library(Quandl)
start_date <- "1990-01-01"; end_training <- "2017-01-01"; 
test <- Quandl("BCB/UDJIAD1",trim_start=start_date)
library(Quandl)
DJIA<-Quandl("BCB/UDJIAD1",trim_start=start_date, trim_end=end_training)
#DJIA <- Quandl("BCB/UDJIAD1")
DJIA <- DJIA[rev(rownames(DJIA)),]
plot(DJIA,type="l")
log_price <- log(DJIA$Value)

# DATA and PARAM
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
data <- list(X=log_price, dt=dt, process=3, scheme=1, jump=11, 
             niter_no_jump=5, niter_jump=64, ghiter_no_jump=16, ghiter_jump=32, 
             line_search = 0)

# Likelihood optimisation
obj <- MakeADFun(data, param)
obj$fn(obj$par)
obj$gr(obj$par)
obj$he(obj$par)
opt <- nlminb(opt$par, obj$fn, obj$gr, control=list(trace=1))
opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he, control=list(trace=1))
obj$fn(opt$par)
res <- sdreport(obj)
res

dyn.unload(dynlib("jd_exact_mspa"))
