############################################################
# dynaical systems approach
############################################################

# probability random graph for 1 given neighbourhood, majority rule
phi <- function(x,n,p=0.23){
	n0 <- round(n/2)
	r <- 0:n
	z <- choose(n,r)*(x^(r))*(1-x)^(n-r)*ifelse(r<=n0,p,1-p)
	phiz <- sum(z)
	return(phiz)
}


# probability random graph for 1 given neighbourhood for Ising probability
psi.rho <- function(x,n,theta=-0.23,theta0=-4,pe=0.2){
	z <- array(0,dim=c(length(x),n))
	nu <- round(n*pe)
	for(i in 0:nu) {
		xn <- 1/(1+exp(-theta0-theta*(i-1)))
		z[,i] <- xn*choose(nu,i-1)*x^(i-1)*(1-x)^(nu-i+1)
	}
	phiz <- apply(z,1,sum)
	return(phiz)
}


############################################################
####  majority process								   


# probability plots for binomial process, majority rule
n <- 20
pe <- 0.3
x <- seq(0,1,len=100)
rhox <- c()
for(i in 1:length(x)) rhox[i] <- phi(x[i],n=round(pe*n),p=0.1)
plot(x,rhox,type='l',bty='n',col='blue',xlab='majority probability at n',ylab='majority probability at n+1',ylim=c(0,1),xlim=c(0,1))
lines(c(0,1),c(0,1),col='gray')
lines(x,rhox,col='blue')
for(i in 1:length(x)) rhox[i] <- phi(x[i],n=round(pe*n),p=0.3)
 for(i in 1:length(x)) rhox[i] <- phi(phi(x[i],n=round(pe*n),p=0.9),n=round(pe*n),p=0.9) # higher order steps
lines(x,rhox,col='red',lty=3)
text(locator(1),expression(theta==0.1),col='blue')
text(locator(1),expression(theta==0.3),col='red')

# probability plots for binomial process, majority rule, higher order steps
n <- 20
pe <- 0.3
x <- seq(0,1,len=100)
rhox <- c()
for(i in 1:length(x)) rhox[i] <- phi(x[i],n=round(pe*n),p=0.9)
plot(x,rhox,type='l',bty='n',col='red',xlab='majority probability at n',ylab='majority probability at n+1/n+2',ylim=c(0,1),xlim=c(0,1),lty=3)
lines(c(0,1),c(0,1),col='gray')
lines(x,rhox,col='red',lty=3)
#for(i in 1:length(x)) rhox[i] <- phi(x[i],n=round(pe*n),p=0.3)
for(i in 1:length(x)) rhox[i] <- phi(phi(x[i],n=round(pe*n),p=0.9),n=round(pe*n),p=0.9) # higher order steps
lines(x,rhox,col='red')
text(locator(1),expression(theta==0.9),col='red')
text(locator(1),expression(f),col='red')
text(locator(1),expression(f^2),col='red')


############################################################
## binomial orbits (majority rule)

# orbit for function psi with Ising model
phi.orbit <- function(x,n,theta,pe=0.2,N=100,M=100,rg=TRUE){
  z <- rep(0,N)
  z[1] <- x
  nu <- round(pe*n)
  for(i in c(2:(N-1))){
	y <- rep(0,nu)
	for(k in 0:nu) {
		xn <- ifelse(k<=nu/2,theta,1-theta)
		if(rg) y[k] <- xn*choose(nu,k)*z[i-1]^(nu-k)*(1-z[i-1])^(k)
		if(!rg) y[k] <- xn
	}
	if(rg) z[i+1] <- sum(y) 
	if(!rg) z[i+1] <- y[nu]
  }
  return(z[c((N-M):N)])
}

# test functions
fo <- phi.orbit(x=runif(1,0,1),theta=0.2,n=20,pe=0.4,rg=TRUE,N=500,M=500)
plot(fo,type='l',col='blue',bty='n',ylim=c(0,1),xlab='time',ylab='Ising probability')

# bifurcation plots
my.r <- seq(0,1, by=0.001)
system.time(Orbit <- sapply(my.r, phi.orbit,  x=runif(1,0,1), n=20,pe=0.3, N=500, M=10))
Orbit <- as.vector(Orbit)
r <- sort(rep(my.r, 10+1))
plot(r,Orbit,pch=16,cex=0.3,col='lightsteelblue4',bty='n',ylim=c(0,1),xlab=expression(theta),ylab='majorityprobability')
lines(c(0.1,0.1),c(0,1),col='blue')
lines(c(0.3,0.3),c(0,1),col='red',lty=3)

############################################################
####  Ising process								       

# plots for Ising process
x <- seq(0,1,len=100)
rhox <- psi.rho(x,n=20,theta=1.2,theta0=-4,pe=0.4)
plot(x,rhox,type='l',bty='n',col='blue',xlab="ising probability at n",ylab="Ising probability at n+1",ylim=c(0,1),xlim=c(0,1))
lines(c(0,1),c(0,1),col='gray')
lines(x,rhox,col='blue')
rhox <- psi.rho(x,n=20,theta=2.7,theta0=-4,pe=0.4)
lines(x,rhox,col='red',lty=3)
text(locator(1),expression(theta==1.2),col='blue')
text(locator(1),expression(theta==2.7),col='red')

# plots for Ising process, higher order steps
x <- seq(0,1,len=100)
rhox <- psi.rho(x,n=20,theta=2.7,theta0=-4,pe=0.4)
plot(x,rhox,type='l',bty='n',col='red',xlab="ising probability at n",ylab="Ising probability at n+1",ylim=c(0,1),xlim=c(0,1),lty=3)
lines(c(0,1),c(0,1),col='gray')
lines(x,rhox,col='red',lty=3)
rhox <- psi.rho(psi.rho(x,n=20,theta=2.7,theta0=-4,pe=0.4),n=20,theta=2.7,theta0=-4,pe=0.4)
lines(x,rhox,col='red')
text(locator(1),expression(f),col='red')
text(locator(1),expression(f^2),col='red')
text(locator(1),expression(theta==2.7),col='red')

############################################################
## Ising orbits

# orbit for function psi with Ising model
psi.orbit <- function(x,n,theta,pe=0.2,N=100,M=100,rg=TRUE,theta0=-4){
  z <- rep(0,N)
  z[1] <- x
  nu <- round(pe*n)
  for(i in c(2:(N-1))){
	y <- rep(0,nu)
	for(k in 1:nu) {
		xn <- 1/(1+exp(-theta*(k-1)-theta0)) 
		if(rg) y[k] <- xn*choose(nu,k-1)*z[i-1]^(k-1)*(1-z[i-1])^(nu-k+1)
		if(!rg) y[k] <- xn
	}
	if(rg) z[i+1] <- sum(y) 
	if(!rg) z[i+1] <- y[nu]
  }
  return(z[c((N-M):N)])
}

# test functions
fo <- psi.orbit(x=runif(1,0,1),theta=1.1,n=20,pe=0.4,rg=TRUE,theta0=-4,N=50000,M=50000)
plot(fo,type='l',col='blue',bty='n',ylim=c(0,1),xlab='time',ylab='Ising probability')

# plot(fo[22001:23000],type='l',col='blue',bty='n',ylim=c(0,1),xlab='time',ylab='Ising probability')

# bifurcation plots
my.r <- seq(0,3, by=0.001)
system.time(Orbit <- sapply(my.r, psi.orbit,  x=runif(1,0,1), n=20, theta0=-4,pe=0.4, N=1000, M=10))
Orbit <- as.vector(Orbit)
r <- sort(rep(my.r, 10+1))
plot(r,Orbit,pch=16,cex=0.3,col='lightsteelblue4',bty='n',ylim=c(0,1),xlab=expression(theta),ylab='Ising probability')
lines(c(1.2,1.2),c(0,1),col='blue')
lines(c(2.7,2.7),c(0,1),col='red',lty=3)


### wave plots of dynamics (a la Palmer)
ktim <- 100
fo <- rep(NA,100)
plot(fo,type='p',col='#0000FF50',bty='n',ylim=c(0,1),xlab='time',ylab='Ising probability',pch=16)
for(i in 1:ktim){
	fo <- psi.orbit(x=runif(1,0,1),theta=0.2,n=20,pe=0.4,rg=TRUE,theta0=-4,N=100,M=100)
	lines(fo,col='#0000FF50') #paletteer_d("colorBlindness::Blue2DarkOrange18Steps")[i])
}


