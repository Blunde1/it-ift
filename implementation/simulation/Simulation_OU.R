#Exact simulation for an OU process
#Author: Berent Å. S. Lunde
#Date: 09.11.2015

OU<-function(dt,kappa,alpha,sigma,x0){
    #dX_t=kappa*(alpha-X_t)*dt + sigma*dW_t
    #X_t|X_0 \sim 
    #constants
    OU.e<-x0*exp(-kappa*dt)+alpha*(1-exp(-kappa*dt))
    OU.v<-sigma^2*(1-exp(-2*kappa*dt))/(2*kappa)
    X<-rnorm(1,mean=OU.e,sd=sqrt(OU.v)) # X = X_{t+dt}
    return(X)
}

OU_process <- function(T,N,kappa,alpha,sigma,x0,seed){
    set.seed(seed)
    x <- c(x0)
    dt<-T/N
    for (i in 2:(N+1)) {
        x[i]  <-  OU(dt,kappa,alpha,sigma,x[i-1])
    }
    return(x);
}
test<-OU_process(10,1000,2,1,1,3,1)
plot(test,type="l")
