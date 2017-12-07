setwd("C:/Users/Berent/Projects/it-ift/implementation")
setwd("~/Projects/UiS-Git/it-ift/implementation")
library(TMB)
compile("jd_ift_quadrature.cpp")
dyn.load(dynlib("jd_ift_quadrature"))

# real data
real_data <- TRUE
if(real_data){
    library(Quandl)
    start_date <- "1950-01-01"; end_training <- "2017-01-01"; 
    test <- Quandl("BCB/UDJIAD1",trim_start=start_date)
    library(Quandl)
    DJIA<-Quandl("BCB/UDJIAD1",trim_start=start_date, trim_end=end_training)
    #DJIA <- Quandl("BCB/UDJIAD1")
    DJIA <- DJIA[rev(rownames(DJIA)),]
    plot(DJIA,type="l")
    log_price <- log(DJIA$Value)
    data <- list(X=log_price, dt=1/250, process=2, scheme=1, jump=0, qiter=100, quadrature=2)
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
    data <- list(X=log(X), dt=1/12, process=2, scheme=1, jump=0, qiter=100, quadrature=2)
    plot(X, type="l", main="Simulated GBM")
}

##### GBM Estimation ####
par_diff <- c(kappa=0.1,sigma=0.2)
par_jump <- c()
param <- list(par = c(par_diff,par_jump))
obj <- MakeADFun(data, param)
opt <- nlminb(obj$par, obj$fn, obj$gr, obj$he, control=list(trace=1))
res <- sdreport(obj)
res

##### MJD Estimation ####
data$process = 3
data$jump = 1
data$qiter = 180
par_jump <- c(30,0,0.05)
param <- list(par = c(par_diff, par_jump))
obj <- MakeADFun(data, param)
opt <- nlminb(obj$par, obj$fn, obj$gr, control=list(trace=1))
res <- sdreport(obj)
res

##### MJD TD Plotting #####
param <- list(par=c(0.4,log(0.3),log(10),-0.01,log(0.05))) # Actual master
# param$par <- c(0.55,log(0.2),log(18),-0.0063,log(0.4))
x0 <- 0; dt <- 1/250
data$X <- rep(x0,2)
data$process = 3
data$jump = 1
data$qiter = 180
data$dt = 1/250

# Matsuda multimodal settings
param <- list(par=c("r"=0.03,"sigma"=log(0.2),"lambda"=log(1),"mu"=-0.5,"nu"=log(0.1)))
data$dt <- 1/4
data$X <- rep(x0<-0, 2)
x_vals <- seq(-1.7,0.5,by=0.01)

x_vals <- seq(-0.12, 0.12, by=0.001)
y_vals <- numeric(length(x_vals))
for(i in 1:length(x_vals)){
    data$X <- c(x0,x_vals[i])
    obj <- MakeADFun(data, param)
    y_vals[i] <- exp(-obj$fn()) # obj returns nll
}
lines(x_vals,y_vals,col="blue")
plot(x_vals, y_vals, type="l") # bimodal works



dyn.unload(dynlib("jd_ift_quadrature"))
