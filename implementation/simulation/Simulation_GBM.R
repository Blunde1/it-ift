#Simulate a GBM

GBM_process <- function(T,N,mu,sigma,x0,seed){
    set.seed(seed)
    dt<-T/N
    W<-rnorm((N-1),(mu-sigma^2/2)*dt,sigma*sqrt(dt))
    lX<-c(log(x0))
    for(i in 2:N){
        lX[i]<-lX[i-1]+W[i-1]
    }
    return(exp(lX))
}
#test
#g<-GBM_process(50,800,0.1,0.2,1,1)
#plot(g,type="l")
#plot(log(g),type="l")
