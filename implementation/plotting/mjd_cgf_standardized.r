# Program for plotting the cgf of a Merton Jump Diffusion

# K

# K'

# K''

# s.hat (spa)

# phi.stand

# plot phi.stand / sqrt(2*pi), versus standard normal

# Plot aeld(exp(-x^2) and phi.stand change of variable)

# K
mjd.cgf <- function(s, r=0.03,sigma=0.2,lambda=1,mu=-0.5,nu=0.1, x0=1, t=1/4){
    k = exp(mu+0.5*nu^2) - 1
    s*x0 + s*t*(r-lambda*k-sigma^2/2) + s^2*sigma^2*t/2 +
        lambda*t*(exp(s*mu+s^2*nu^2/2)-1)
}
# K'
mjd.cgf.1 <- function(s, r=0.03,sigma=0.2,lambda=1,mu=-0.5,nu=0.1, x0=1, t=1/4){
    k = exp(mu+0.5*nu^2) - 1
    x0+t*(r-lambda*k-sigma^2/2) + s*sigma^2*t + 
        (mu+nu^2*s) * lambda*t*exp(s*mu+s^2*nu^2/2)
}
# K''
mjd.cgf.2 <- function(s, r=0.03,sigma=0.2,lambda=1,mu=-0.5,nu=0.1, x0=1, t=1/4){
    k = exp(mu+0.5*nu^2) - 1
    sigma^2*t + lambda*t*exp(s*mu+s^2*nu^2/2) * (nu^2+(mu+nu^2*s)^2)
}
# phi.stand
s.hat.fn <- function(s.start, x, cgf=mjd.cgf, cgf.1=mjd.cgf.1){
    nlminb(start=s.start, 
           objective=function(s,x){cgf(s)-s*x}, 
           gradient=function(s,x){cgf.1(s)-x}, 
           x=x)$par
}
# phi.stand
phi.stand <- function(s, x, cgf=mjd.cgf, cgf.1=mjd.cgf.1, cgf.2=mjd.cgf.2){
    s.hat <- s.hat.fn(0, x, cgf, cgf.1)
    exp(
        -cgf(s.hat) - 1i*s*x/sqrt(cgf.2(s.hat)) +
            cgf(
                1i*s/sqrt(cgf.2(s.hat)) + s.hat
            )
    )
}
phi.stand(0, 10)

cgf.stand <- function(s, x, cgf=mjd.cgf, cgf.1=mjd.cgf.1, cgf.2=mjd.cgf.2){
    s.hat <- s.hat.fn(0, x, cgf, cgf.1)
    (
        -cgf(s.hat) - s*x/sqrt(cgf.2(s.hat)) +
            cgf(
                s/sqrt(cgf.2(s.hat)) + s.hat
            )
    )
}
cgf.stand(0,10)
check <- nlminb(10, cgf.stand, x=100)$par
s <- seq(-6,6, length.out = 100)
y.cgf <- sapply(s, cgf.stand, x=12)
plot(s,y.cgf, type="l")

# plot phi.stand / sqrt(2*pi), versus standard normal
u <- seq(-6,6,0.01)
y <- sapply(u, function(u){Re(phi.stand(u, x=10))/(sqrt(2*pi))})
plot(u,y, type="l")
lines(u, dnorm(u,0,1),col="red")
legend("topright", legend=c("mjd standardized", "standard normal"), col=c("black", "red"), lty=c(1,1))

# as a function of x
p.0 <- function(x){
    u <- seq(-6,6,0.05)
    y <- sapply(u, function(u){Re(phi.stand(u, x))/(sqrt(2*pi))})
    sum(y)*(u[2]-u[1])
}
mjd.e <- mjd.cgf.1(0)
mjd.v <- mjd.cgf.2(0)
x <- seq(mjd.e - 2.5*sqrt(mjd.v), mjd.e+2.5*sqrt(mjd.v), length.out = 100)
p.0.x <- sapply(x, p.0)
plot(x, p.0.x,type="l")


# Plot aeld(exp(-x^2) and phi.stand change of variable)



# Approximate int phi.stand / sqrt(2*pi) with spa
s.stand.hat <- function(start, cgf, x){
    nlminb(start, cgf, x=x)$par
}

cgf.stand.1 <- function(s, x, cgf=mjd.cgf, cgf.1=mjd.cgf.1, cgf.2=mjd.cgf.2){
    s.hat <- s.hat.fn(0, x, cgf, cgf.1)
    (
         - x/sqrt(cgf.2(s.hat)) +
            cgf.1(
                s/sqrt(cgf.2(s.hat)) + s.hat
            ) / sqrt(cgf.2(s.hat))
    )
}
cgf.stand.2 <- function(s, x, cgf=mjd.cgf, cgf.1=mjd.cgf.1, cgf.2=mjd.cgf.2){
    s.hat <- s.hat.fn(0, x, cgf, cgf.1)
    (
            cgf.2(
                s/sqrt(cgf.2(s.hat)) + s.hat
            ) / cgf.2(s.hat)
    )
}
s.stand.hat.0 <- s.stand.hat(0, cgf.stand, x=12)
cgf.stand.1(s.stand.hat.0, x=12)
cgf.stand.2(s.stand.hat.0, x=12)

spa <- function(sp, x, cgf=mjd.cgf, cgf.2=mjd.cgf.2){
    exp(cgf(sp)-sp*x) / sqrt(2*pi*cgf.2(sp))
}
doublespa <- function(spa1, sp.stand, x){
    spa1 * exp(cgf.stand(sp.stand, x)) / sqrt(cgf.stand.2(sp.stand, x))
}
x = 0
sp <- s.hat.fn(0,x)
spa1 <- spa(sp, x)
sp.stand <- s.stand.hat(0, cgf.stand, x=x)
doublespa(spa1, sp.stand, x)

x <- seq(-0.75,0.75,by=0.01)
spa1 <- sapply(x, function(x){
    sp <- s.hat.fn(0,x)
    spa(sp, x)
})
spa2 <- sapply(x, function(x){
    sp <- s.hat.fn(0,x)
    spa1 <- spa(sp, x)
    sp.stand <- s.stand.hat(0, cgf.stand, x=x)
    doublespa(spa1, sp.stand, x)
})
p.0.x <- sapply(x, p.0)

plot(x,spa1, ylim=c(0,4), type="l")
plot(x, spa2)
plot(x,spa1*p.0.x, type="l")

# Test for finding spa
x <- 10
s <- seq(-5,5,by=0.0001)
y.inner <- sapply(s, function(s){mjd.cgf(s)-s*x})
y.1 <- sapply(s, function(s){mjd.cgf.1(s)-x})
plot(s,y.inner)
abline(v=s[which.min(y.inner)])
abline(v=s[which.min(abs(y.1))], col="red")
abline(v=s.hat(s.start = 0, x=x)$par, col="blue")

