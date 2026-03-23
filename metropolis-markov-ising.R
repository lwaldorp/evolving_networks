#####################################################
# Posterior sampling mcmc using Metropolis algorithm
#####################################################

# prior: Gamma(3,1)
# likelihood: prod(sum_r(1/(1+exp(-theta0-theta*r)))*choose(nu,r)*x_n^r*(1-x_n)^(nu-r))
# posterior: likelihood * prior


# likelihood
# probability random graph for 1 given neighbourhood for Ising probability
psi.rho <- function(x,n,theta=-0.23,rg=TRUE,pe=0.2,theta0=3){
	z <- array(0,dim=c(length(x),n))
	for(i in 1:n) {
		xn <- 1/(1+exp(-theta0-theta*(i-1))) #+runif(1,0,0.2)
		if(rg) 	 z[,i] <- xn*choose(n,i-1)*x^(i-1)*(1-x)^(n-i+1)
		if(!rg)  z[,i] <- xn
	}
	if(rg) phiz <- apply(z,1,sum)
	if(!rg) phiz <- z[,n] 
	phiz[phiz==0] <- 0.0001
	return(phiz)
}



# metropolis algorithm for both variable, theta0 and theta
metropolis <- function(theta,n=20,x,theta0=-2,pe=0.4,alpha=2,alpha0=3,rate=1,scale=1) {
  theta.new <- runif(1,theta-alpha,theta+alpha)
  theta0.new <- runif(1,theta0-alpha0,theta0+alpha0)
  post.new <- sum(log(psi.rho(x,n=n,theta=theta.new,theta0=theta0.new,pe=pe))) + log(dgamma(theta.new,rate,scale)*dgamma(theta0.new,rate,scale)+0.0001) #+ log(dunif(theta,theta-alpha,theta+alpha)*dunif(theta0,theta0-alpha0,theta0+alpha0)+0.0001)# log(likelihood) + log(prior)
  post.cur <- sum(log(psi.rho(x,n=n,theta=theta,theta0=theta0,pe=pe))) + log(dgamma(theta,rate,scale)*dgamma(theta0,rate,scale)+0.0001) #+ log(dunif(theta.new,theta-alpha,theta+alpha)*dunif(theta0.new-alpha0,theta0.new+alpha0)+0.0001)
  lrat <- post.new - post.cur
  if(log(runif(1)) < lrat){
  	theta <- theta.new
  	theta0 <- theta0.new
  	}
  return(c(theta=theta,theta0=theta0))
}


# settings
set.seed(1492)
nrun <- 10000
x <- numeric()
y <- numeric()
theta.0 <- runif(1,0,2)
theta0.0 <-  runif(1,0,2)
theta.cur <- metropolis(theta=theta.0,theta0=theta0.0,x=rho,alpha=1,alpha0=1)
for(j in 1:nrun) {
  theta.cur <- metropolis(theta=theta.cur[1],theta0=theta.cur[2], x=rho,alpha=1,alpha0=1,mu1=1,sigma1=1,mu=-2,sigma=1,pe=0.4)
  x[j] <- theta.cur[1]
  y[j] <- theta.cur[2]
}

# select after burn-in
burnin <- 5000
x <- x[burnin:length(x)]
y <- y[burnin:length(y)]

# checks for theta
plot(x,col="#0000FF70",bty='n',xlab='iteration',ylab=expression(theta),type='l')
hist(x,border='white',col='#0000FF40',bty='n',xlab=expression(theta),ylab='density',main='',probability=TRUE)
lines(density(x),bty='n',col='#0000FF70')
d <- density(x)
mean(x); d$x[which.max(d$y)]

# checks for theta0
plot(y,col="#0000FF70",bty='n',xlab='iteration',ylab=expression(theta[0]),type='l')
hist(y,border='white',col='#0000FF40',bty='n',xlab=expression(theta[0]),ylab='density',main='',probability=TRUE)
lines(density(y),bty='n',col='#0000FF70')
d <- density(y)
mean(y); d$x[which.max(d$y)]

# plot confidence interval in density
# library(bayesplot)
Z <- cbind(theta=x,theta0=y)
mcmc_areas(Z,pars="theta",prob=0.95)


# Groot data
# library(mgm)
groot <- symptom_data #in package mgm
mood <- apply(groot$data[,c(1:12)],1,sum) # mood
symp <- apply(groot$data[,c(13:16)],1,sum) # symptoms, agitate, worry etc
se <- apply(groot$data[,c(17:20)],1,sum) # self esteem
soc <- apply(groot$data[,c(21:28)],1,sum) # social activity
phys <- apply(groot$data[,c(29:33)],1,sum) # tired, hyngry etc
act <- apply(groot$data[,c(34:38)],1,sum) # tired, hyngry etc
event <- apply(groot$data[,c(40:47)],1,sum) # pleasure, etc 

groot.score <- cbind(mood,symp,se,phys)
colnames(groot.score) <- labs <- c("mood","symp","se","phys")
groot.sum <- apply(groot.score,1,sum)
xg <- groot.sum/max(groot.sum)
plot(xg,bty='n',col='#0000FF60',type='l',xlab='time',ylab='proportion',ylim=c(0.5,1))
hist(xg,bty='n',border='#FFFFFF',col='#0000FF70',main='',prob=TRUE,xlab='proportion',ylab='density')
lines(density(xg),col='#0000FF80')

# Bayes estimation 
nrun <- 10000
x <- numeric()
y <- numeric()
theta.0 <- 0.5
theta0.0 <- -0.5
theta.cur <- metropolis(theta=theta.cur[1],theta0=theta.cur[2], x=xg,alpha=1,alpha0=1,mu1=0.7,sigma1=1,mu=-2,sigma=1,pe=0.4)
# metropolis(theta=theta.0,theta0=theta0.0,x=xg,alpha=2,alpha0=5,pe=0.3)
for(j in 1:nrun) {
  theta.cur <- metropolis(theta=theta.cur[1],theta0=theta.cur[2], x=rho,alpha=1,alpha0=1)
  x[j] <- theta.cur[1]
  y[j] <- theta.cur[2]
}

# select after burn-in
burnin <- 5000
x <- x[burnin:length(x)]
y <- y[burnin:length(y)]

# checks for theta
plot(x,col="#0000FF70",bty='n',xlab='iteration',ylab=expression(theta),type='l')
hist(x,border='white',col='#0000FF40',bty='n',xlab=expression(theta),ylab='density',main='',probability=TRUE)
lines(density(x),bty='n',col='#0000FF70')
d <- density(x)
mean(x); d$x[which.max(d$y)]
# 95% credible interval, complete time series
# sum(d$y[d$x<=0.34])/sum(d$y) # 0.02507812
# sum(d$y[d$x>=3.2])/sum(d$y) # 0.02493615
# mode = 1.253


# checks for theta0
plot(y,col="#0000FF70",bty='n',xlab='iteration',ylab=expression(theta[0]),type='l')
hist(y,border='white',col='#0000FF40',bty='n',xlab=expression(theta[0]),ylab='density',main='',probability=TRUE)
lines(density(y),bty='n',col='#0000FF70')
d <- density(y)
mean(y); d$x[which.max(d$y)]
# 95% credible interval, complete time series
# sum(d$y[d$x<=-3.35])/sum(d$y) # 0.02541239
# sum(d$y[d$x>=4.4])/sum(d$y) # 0.02545926
# mode = 0.4435

