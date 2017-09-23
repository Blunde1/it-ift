# normal ift R
# Berent Lunde
# 21.09.2017

setwd("~/Projects/UiS-Git/it-ift/implementation")
library(TMB)
compile("normal_ift.cpp")
dyn.load(dynlib("normal_ift"))

data <- list(x=1, gqiter=150) # problems in tails, solved if gqiter is high
param <- list(mu=1,sigma=2)
obj <- MakeADFun(data,param)
obj$report()
obj$fn(obj$par)

f <- x <- with(param,seq(mu-6*sigma,mu+6*sigma,length.out = 100))
for(i in 1:length(x)){
    data$x <- x[i]
    obj <- MakeADFun(data,param)
    f[i] <- obj$fn(obj$par)
}
with(param,plot(x,dnorm(x,mu,sigma),type="l"))
lines(x,f, lty=3,col="red")
exact <- with(param, dnorm(x,mu,sigma))
plot(x,abs(log(exact)-log(f)),type="l")

dyn.unload(dynlib("normal_ift"))
