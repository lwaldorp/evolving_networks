############################################################
# Markov chain approach
############################################################

# probability random graph for 1 given neighbourhood for majority probability
phi <- function(x,n,p=0.23){
	n0 <- round(n/2)
	r <- 0:n0
	z <- choose(n,r)*(x^(r))*(1-x)^(n-r)*ifelse(r<n0,p,1-p)
	phiz <- sum(z)
	return(phiz)
}


# probability random graph for 1 given neighbourhood for Ising probability
psi.rho <- function(x,n,theta=-0.23,rg=TRUE,pe=0.2,theta0=3){
	nu <- round(n*pe)
	z <- array(0,dim=c(length(x),nu))
	for(i in 0:nu) {
		xn <- 1/(1+exp(-theta0-theta*(i-1))) 
		if(rg) 	 z[,i] <- xn*choose(nu,i-1)*x^(i-1)*(1-x)^(nu-i+1)
		if(!rg)  z[,i] <- xn
	}
	if(rg) phiz <- apply(z,1,sum)
	if(!rg) phiz <- z[,n] 
	return(phiz)
}


############################################################
####  majority process								    ####
############################################################
set.seed(1492)
nrun <- 10
rho_nrun <- array(NA,dim=c(nsim,nrun))
for(k in 1:nrun){
n <- 20
p <- 0.2
pe <- 0.4
nsim <- 100
nrho <- c(1:nsim)
nrho[1] <- runif(1,0,1) #round(p*n,0)
for(i in 2:nsim) nrho[i] <- rbinom(1,n,phi(nrho[i-1]/n,round(pe*n,0),p)) #rbinom(1,n,phi(nrho[i-1]/n,n,p))

# plot(rho,type='l',col='blue',bty='n',xlab='time',ylab="proportion on",axes=TRUE,ylim=c(0,1))	

rho <- nrho/n
rho_nrun[,k] <- rho
}

rho_mean <- apply(rho_nrun,1,mean)
rho_sd   <- apply(rho_nrun,1,sd)
rho_lb   <- rho_mean-2*rho_sd
rho_ub   <- rho_mean+2*rho_sd

plot(rho_mean, type="l", lty=1, bty="n", col='blue', xlab="time", ylab="proportion on nodes",ylim=c(0,1))
polygon(c(1:nsim, rev(1:nsim)), c(rho_ub, rev(rho_lb)), col = "#0000FF30", border=0)
lines(rho, type="l", col="red", lty=3)

############################################################
####  Ising process             					    ####
############################################################
set.seed(1492)
nrun <- 10
rho_nrun <- array(NA,dim=c(nsim,nrun))
n <- 20
theta0 <- 0
theta <- 1.2
pe <- 0.4
nsim <- 100

for(k in 1:nrun){
nrho <- c(1:nsim)
nrho[1] <- runif(1,0,1)
for(i in 2:nsim) nrho[i] <- rbinom(1,n,psi.rho(nrho[i-1]/n,n=round(pe*n,0),theta=theta, theta0=theta0,pe=pe)) #rbinom(1,n,phi(nrho[i-1]/n,n,p))

# plot(rho,type='l',col='blue',bty='n',xlab='time',ylab="proportion on",axes=TRUE,ylim=c(0,1))	

rho <- nrho/n
rho_nrun[,k] <- rho
}

rho_mean <- apply(rho_nrun,1,mean)
rho_sd   <- apply(rho_nrun,1,sd)
rho_lb   <- rho_mean-2*rho_sd
rho_ub   <- rho_mean+2*rho_sd

plot(rho_mean, type="l", lty=1, bty="n", col='blue', xlab="time", ylab="proportion on nodes",ylim=c(0,1))
polygon(c(1:nsim, rev(1:nsim)), c(rho_ub, rev(rho_lb)), col = "#0000FF30", border=0)
lines(rho, type="l", col="red", lty=3)
