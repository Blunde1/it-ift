# Program for inspacting the cgf of a Merton Jump Diffusion

param <- c("r"=0.03,"sigma"=0.2,"lambda"=1,"mu"=-0.5,"nu"=0.1, "x0"=0, "t"=1/4)

# K
mjd.cgf <- function(s, r=param["r"],sigma=param["sigma"],lambda=param["lambda"],
                    mu=param["mu"],nu=param["nu"], x0=param["x0"], t=param["t"]){
    k = exp(mu+0.5*nu^2) - 1
    s*x0 + s*t*(r-lambda*k-sigma^2/2) + s^2*sigma^2*t/2 +
        lambda*t*(exp(s*mu+s^2*nu^2/2)-1)
}
# K'
mjd.cgf.1 <- function(s, r=param["r"],sigma=param["sigma"],lambda=param["lambda"],
                      mu=param["mu"],nu=param["nu"], x0=param["x0"], t=param["t"]){
    k = exp(mu+0.5*nu^2) - 1
    x0+t*(r-lambda*k-sigma^2/2) + s*sigma^2*t + 
        (mu+nu^2*s) * lambda*t*exp(s*mu+s^2*nu^2/2)
}
# K''
mjd.cgf.2 <- function(s, r=param["r"],sigma=param["sigma"],lambda=param["lambda"],
                      mu=param["mu"],nu=param["nu"], x0=param["x0"], t=param["t"]){
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
y.cgf <- sapply(s, cgf.stand, x=10)
plot(s,y.cgf, type="l")

# plot phi.stand / sqrt(2*pi), versus standard normal
u <- seq(-10,10,0.01)
y <- sapply(u, function(u){Re(phi.stand(u, x=0))/(sqrt(2*pi))})
plot(u,y, type="l", ylim=c(min(y), max(y)), col="1", lwd=2, lty=2, xlab="s", ylab=expression(paste(varphi[hat(x)(tilde(tau))](s))))
lines(u, dnorm(u,0,1),col="red")
lines(u, abs(dnorm(u,0,1, log = TRUE) - log(y)), col="blue")
legend("topright", legend=c("mjd standardized", "standard normal"), col=c("black", "red"), lty=c(1,1))

##### PLOTTING W.R.T. S
pdf("mjd_stand_s.pdf", width=7, height=4+1/3)

par(mar=c(5,6, 3, 3))
plot(u, y, xaxt="n", yaxt="n", ylim=c(min(y), max(y)), 
     , main="", type="l", col="1", lwd=2, lty=2,
     ylab=expression(paste(Re,"[", varphi[hat(x)(tilde(tau))](s),"]")),
     xlab = expression(paste(s)))
lines(u, dnorm(u,0,1), col=2, lwd=2, lty=2)
axis(2, ylim=c(0, max(y)), col="black", lwd=1, las=1)
axis(1,ylab=expression(s))
legend("topright", legend=c("MJD standardised", "Standard normal"), col=c("black", "red"), lty=c(2,2))

dev.off()



##### PLOTTING W.R.T. X


# as a function of x
p.0 <- function(x){
    u <- seq(-10,10,0.02)
    y <- sapply(u, function(u){Re(phi.stand(u, x))/(2*pi)})
    sum(y)*(u[2]-u[1])
}
mjd.e <- mjd.cgf.1(0)
mjd.v <- mjd.cgf.2(0)
x <- seq(mjd.e - 9*sqrt(mjd.v), mjd.e+2*sqrt(mjd.v), length.out = 200)
p.0.x <- sapply(x, p.0)
plot(x, p.0.x,type="l")
spa <- function(sp, x, cgf=mjd.cgf, cgf.2=mjd.cgf.2){
    exp(cgf(sp)-sp*x) / sqrt(2*pi*cgf.2(sp))
}
spa1 <- sapply(x, function(x){
    sp <- s.hat.fn(0,x)
    spa(sp, x)
})

pdf("mjd_prob_x.pdf", width=7, height=4+1/3)
par(mar=c(5,6, 3, 6))
plot(x, spa1, xaxt="n", yaxt="n", ylim=c(0, max(spa1*p.0.x*sqrt(2*pi))), 
     , main="", type="l", col="1", lwd=2, lty=2,
     ylab=expression(paste("Spa(f;x)")),
     xlab = expression(paste(x)))
lines(x, spa1*p.0.x*sqrt(2*pi), col=2, lwd=2, lty=2)
axis(2, ylim=c(0, max(spa1*p.0.x*sqrt(2*pi))), col="black", lwd=1, las=1)
axis(1)
par(new=T)
plot(x, p.0.x, axes=F, xlab="", ylab="", main="", type="l", col=4, lty=2, lwd=1.5)
axis(4, ylim=c(min(p.0.x), max(p.0.x)))
mtext(4, text=expression(paste(
    p[bar(X)(hat(tau))](0)
)), line=3)

legend("topleft", legend=c("Orindary spa", "Exact", "Probability in zero"), col=c("black", "red", 4), 
       lty=c(2,2, 2), lwd=c(2,2, 1.5))
dev.off()

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
lines(x, spa2, col="red")
lines(x,spa1*p.0.x, type="l", col="blue")

# Test for finding spa
x <- 10
s <- seq(-5,5,by=0.0001)
y.inner <- sapply(s, function(s){mjd.cgf(s)-s*x})
y.1 <- sapply(s, function(s){mjd.cgf.1(s)-x})
plot(s,y.inner)
abline(v=s[which.min(y.inner)])
abline(v=s[which.min(abs(y.1))], col="red")
abline(v=s.hat(s.start = 0, x=x)$par, col="blue")

