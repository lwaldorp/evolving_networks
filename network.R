############################################################
# network of ising/majority nodes 
############################################################

library(igraph)

# parameters to set
prob_init <- 0.2 # initial probability of 0/1 nodes
theta0 <- -4 # offset parameter Ising model
theta <- 1.2 # parameter for probaility models
N <- 100 # number of discrete time steps to evolve
m <- 20 # size of the network
pe <- 0.4 # probability of an edge in the random graph

# majority rule
maj_n <- function(x, theta){
	nsize <- length(x)
	r <- sum(x)
	maj <- ifelse(r <= nsize/2, theta, 1 - theta)
	return(maj)
}


# Ising rule
ising_n <- function(x, theta, theta0=-4){
	r <- sum(x)
	ising <- 1/(1 + exp(-theta0 - r*theta))
	return(ising)
}

# g_theta determines 0/1 from probability
g_theta <- function(prob_theta){
	U <- runif(1,0,1)
	g <- ifelse(prob_theta >= U, 1, 0)
	return(g)
}

# create random graph
graph <- erdos.renyi.game(m, p.or.m=pe, type="gnp", directed=FALSE, loops=FALSE)
plot(graph,  layout=layout_on_grid, vertex.label=NA, vertex.size=30, vertex.frame.color="#0000FF50", vertex.label.dist=0.5, vertex.color="#0000FF50", vertex.shape="circle", edge.width=2, rescale=TRUE, asp=1)

# determine set of neighbourhoods in graph 
neighbourhoods <- ego(graph, order=1)

## single run

# initial state, random
x_init <- sample(0:1, size=m, replace=TRUE, prob=c(prob_init, 1-prob_init))

# evolve network for N steps
theta <- rnorm(m,mean=0,sd=1)
rule <- "ising"
x <- x_init # initial value at n=0
x_N <- array(NA,dim=c(N,m))
for(i in 1:N){
	x_plus1 <- c()
	for(j in 1:m){ # no difference if ordered randomly
		if(rule=="ising") x_plus1[j] <- g_theta(ising_n(x[neighbourhoods[[j]]], theta[j], theta0))
		if(rule=="maj")   x_plus1[j] <- g_theta(maj_n(x[neighbourhoods[[j]]], theta))
	}
	x_N[i,] <- x_plus1 
}

x_N <- rbind(x_init,x_N)

evo_prop <- apply(x_N, 1, sum)/m

# plot evolutions
plot(0:N, evo_prop, type="l", bty="n", col='blue', xlab="time", ylab="proportion on nodes",ylim=c(0,1))



## multiple runs
# parameters to set
prob_init <- runif(1,0,1) #0.2 # initial probability of 0/1 nodes
theta0 <- -4 # offset parameter Ising model
theta <- 1.2 # parameter for probaility models
N <- 100 # number of discrete time steps to evolve
m <- 20 # size of the network
pe <- 0.4 # probability of an edge in the random graph

nrun <- 100
evo_prop_nrun <- array(NA,dim=c(N+1,nrun))

for(k in 1:nrun){
# initial state, random
	prob_init <- 0.3 #runif(1,0,1) #0.2 # initial probability of 0/1 nodes
	x_init <- sample(0:1, size=m, replace=TRUE, prob=c(prob_init, 1-prob_init))

# evolve network for N steps
	 theta <- rnorm(m,mean=1.2,sd=1) # draw theta from normal distribution
	# theta <- 1.2 # theta fixed
	# theta <- rgamma(m,shape=1,rate=1)
	rule <- "ising"
	x <- x_init # initial value at n=0
	x_N <- array(NA,dim=c(N,m))
	for(i in 1:N){
	x_plus1 <- c()
	for(j in 1:m){
		if(rule=="ising") x_plus1[j] <- g_theta(ising_n(x[neighbourhoods[[j]]], theta[j], theta0))
		if(rule=="maj")   x_plus1[j] <- g_theta(maj_n(x[neighbourhoods[[j]]], theta))
		}
		x_N[i,] <- x_plus1 
	}

	x_N <- rbind(x_init,x_N)

	evo_prop <- apply(x_N, 1, sum)/m

	evo_prop_nrun[,k] <- evo_prop

	# plot evolutions
	if(k==1) plot(0:N, evo_prop, type="l", lty=1, bty="n",col='#0000FF50',xlab="time",ylab="proportion on nodes",ylim=c(0,1))
else lines(0:N, evo_prop, type="l", col="#0000FF50", lty=1)

}

evo_mean <- apply(evo_prop_nrun,1,mean)
evo_sd   <- apply(evo_prop_nrun,1,sd)
evo_lb   <- evo_mean-2*evo_sd
evo_ub   <- evo_mean+2*evo_sd

plot(0:N, evo_mean, type="l", lty=1, bty="n", col='blue', xlab="time", ylab="proportion on nodes",ylim=c(0,1))
# lines(0:N, evo_mean-2*evo_sd, col="gray")
# lines(0:N, evo_mean+2*evo_sd, col="gray")
polygon(c(0:N, rev(0:N)), c(evo_ub, rev(evo_lb)), col = "#0000FF30", border=0)
lines(0:N, evo_prop, type="l", col="red", lty=3)

text(locator(1),"Gamma(1,1)")
text(locator(1),"N(1.2,1)")