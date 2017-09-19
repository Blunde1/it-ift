CIR<-function(dt,kappa,alpha,sigma,x0){
    #dX_t=kappa*(alpha-X_t)*dt + sigma*sqrt(X_t)*dW_t
    #X_t|X_0 \sim 
    #constants
    c<-2*kappa/(sigma^2*(1-exp(-kappa*dt)))
    Y<-rchisq(n=1,df=4*kappa*alpha/sigma^2,ncp=2*c*x0*exp(-kappa*dt))
    X<-Y/(2*c) # X = X_{t+dt}
    return(X)
}

CIR_process <- function(T,N,kappa,alpha,sigma,x0,seed){
    set.seed(seed)
    x <- c(x0)
    dt<-T/N
    for (i in 2:(N+1)) {
        x[i]  <-  CIR(dt,kappa,alpha,sigma,x[i-1])
    }
    return(x);
}
#test<-CIR_process(10,1000,2,1,1,3,2)
#plot(test,type="l")
