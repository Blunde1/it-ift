#Exact simulation for an MJD process
#Author: Berent Å. S. Lunde

mjd_process<-function(N,T,r,sigma,lambda,mu,nu,x0,seed){
    set.seed(seed)
    time <- seq(0,T,length.out=N)
    lX<-c(log(x0))
    dt <- T/N
    k <- exp(mu+0.5*nu^2) - 1
    W<-rnorm((N-1),0,sqrt(dt))
    jtimes<-c(0)
    
    for(i in 2:N){
        N <- rpois(1,dt*lambda)
        if(N != 0){ jtimes[length(jtimes)+1] <- i  }
        jump <- sum(rnorm(N,mu,nu))
        lX[i] = lX[i-1] + (r-0.5*sigma^2-lambda*k)*dt +
            sigma*W[i-1] +  jump
    }
    
    res <- list()
    res$x <- exp(lX)
    res$time <- time
    res$jtimes <- jtimes
    return(res)
#    return(exp(lX))
}

test <- mjd_process(3*1000,3,0.4,0.3,30,-0.01,0.05,100,1)
attach(test)

plot(time[1:(jtimes[2]-1)],x[1:(jtimes[2]-1)],xlim=range(time),ylim=range(x),
     type="l", las=1, lwd=1,lty=1,
     xlab="Time",ylab="Value")
for(i in 2:(length(jtimes)-1)){
    lines(time[jtimes[i]:(jtimes[i+1]-1)], x[jtimes[i]:(jtimes[i+1]-1)] ,lwd=1)
}
lines(time[jtimes[length(jtimes)]:length(time)],x[jtimes[length(jtimes)]:length(time)],lwd=1)

obs <- x[seq(1,length(x),length.out=36)]
obstime <- seq(0,3,length.out=36)
points(obstime,obs,pch=4,col="red",lwd=2)
legend("topleft",c("Underlying jump-diffusion process","Observations"),
       col=c("black","red"),lwd=c(1,3),pch=c(NA,4),lty=c(1,NA))

#plot(log(test),type="l")
#c(0.4,0.3,30,-0.01,0.05)