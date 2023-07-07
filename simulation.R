################
### Packages ###
################

require(igraph)
require(Pareto)


##############
#### Setup ###
##############

## Weak dependence if applied: Variogram matrix is multiplied by 4.
## Thus, random fluctuations are more due to the log-normal part
## and less due to the Pareto variable.
weak_dependence = FALSE

## Two-sided jumps if applied: The Lévy measure has mass on all orthants
## of R^d. Otherwise each components has only positive jumps.
two_sided = TRUE

## Dimension of the process
d <- 12

## Stability index
alpha <- 3/2

## Marginal masses m_1,...,m_d (up to a constant)
m <- rep(1,d)


#########################
### Defining the tree ###
#########################

## Matrix of edges in the tree
tree_mat <- matrix(c(
  1,2,
  1,3,
  1,4,
  4,5,
  1,6,
  6,7,
  6,8,
  1,9,
  9,10,
  10,11,
  10,12
),ncol=2,byrow=T)

## Creating and plotting the tree
Tree <- graph.data.frame(tree_mat[,1:2],directed=F)
plot(Tree,asp=0)


##########################
### Hüsler-Reiss setup ###
##########################

## Square root of the matrix S
A <- matrix(c( 
  1,0,0,0,0,0,0,0,0,0,0,0,
  1,1,0,0,0,0,0,0,0,0,0,0,
  1,0,1,0,0,0,0,0,0,0,0,0,
  1,0,0,1,0,0,0,0,0,0,0,0,
  1,0,0,1,1,0,0,0,0,0,0,0,
  1,0,0,0,0,1,0,0,0,0,0,0,
  1,0,0,0,0,1,1,0,0,0,0,0,
  1,0,0,0,0,1,0,1,0,0,0,0,
  1,0,0,0,0,0,0,0,1,0,0,0,
  1,0,0,0,0,0,0,0,1,1,0,0,
  1,0,0,0,0,0,0,0,1,1,1,0,
  1,0,0,0,0,0,0,0,1,1,0,1
),nrow=d,byrow=T)

## Adjusts A if week dependence is selected
if (weak_dependence) {
  A <- 2 * A
}

## Covariance matrix
S <- A%*%t(A)

## Calculating variogram matrix G corresponding to S
G <- cbind(rep(1,d))%*%diag(S)+cbind(diag(S))%*%rep(1,d)-2*S


#############################################
### Parameters for approximate simulation ###
#############################################

## epsilon threshold to be used for approximate simulation
epsilon <- 10^(-3)

## Average number of jumps larger than epsilon
avg_no_jumps <- 10000


########################################
### Estimation of certain parameters ###
########################################

## Function for simulating n jumps
sim_jump <- function(n) {
  res <- matrix(nrow=0,ncol=d)
  while (nrow(res)<n) {
    P <- rPareto(n,epsilon,alpha)
    U <- A%*%matrix(rnorm(n*d),nrow=d)
    direction <- sample(1:d,n,replace=T,prob=m)
    temp <- sapply(1:n,function(k) alpha*log(P[k])+log(m)-U[direction[k],k]-log(m[direction[k]])-1/2*G[,direction[k]])
    X <- exp(1/alpha*(U+temp))
    unif <- runif(n)
    p <- (colSums(X>=epsilon))^(-1) #Acceptance probabilities
    res <- rbind(res,t(X[,unif<p]))
  }
  if (two_sided) {
    signs <- matrix(sample(c(-1,1),n*d,replace=T),ncol=d)
    res = signs*res[1:n,]
  }
  res[1:n,]
}

## Estimation of constant in density of marginals
R1 <- 10^6
X <- sim_jump(R1)
C <- avg_no_jumps*mean(abs(X[,1])>epsilon)*epsilon^alpha*alpha

rm(R1, X)

## Estimation of covariance matrix for small-jump approximation
expectation <- function(U,x) {
  n <- ncol(U)
  matrix_sum <- matrix(0,nrow=d,ncol=d)
  temp <- sapply(1:n,function(k) log(m)-U[1,k]-log(m[1])-1/2*G[,1])
  Y <- exp(1/alpha*(U+temp))
  sapply(1:n,function(k) {matrix_sum <<- matrix_sum + Y[,k]%*%t(Y[,k])*(x*max(Y[,k])<=1)})
  matrix_sum/n
}

R2 <- 10^5
U <- A%*%matrix(rnorm(R2*d),nrow=d)

x_seq <- seq(from=10^(-3),to=1,by=10^(-3))
riemann_sum <- matrix(0,nrow=d,ncol=d)
sapply(x_seq,function(x) {riemann_sum <<- riemann_sum + x^(1-alpha)*expectation(U,x)})
# Covariance matrix for small-jump approximation
Sigma <- C*riemann_sum/length(x_seq)

rm(R2, U)

## Estimation of drift

if (two_sided) {
  g <- rep(0,d)
} else {
  ## alpha>1:
  R1 <- 10^7
  X <- sim_jump(R1)
  g <- -colMeans(X)*avg_no_jumps
  rm(R1, X)
  ## For alpha<1 a Riemann sum approximation is needed (like for the
  ## covariance matrix Sigma)
}


##################
### Simulation ###
##################

## Simulation of the compound Poisson process.
## n is number of subintervals, rate is that of Poisson process (i.e. avg. no. of jumps)
sim_CPP <- function(n,rate) { 
  ## Number of jumps in subintervals
  N <- rpois(n,rate/n) 
  ## Total number of jumps
  M <- sum(N)
  J <- sim_jump(M)
  cJ <- apply(J,2,cumsum)
  cN <- cumsum(N)
  X <- cJ[cN,]
  #Making sure X has exactly n rows
  X <- rbind(matrix(rep(0,d*(n-nrow(X))),ncol=d),X)
  X
}

## Approximate simulation of the Lévy process
##n is number of subintervals
sim_X <- function(n) {
  ## Linear drift
  drift <- (1:n)%*%rbind(g)/n
  ## Brownian motion
  U <- epsilon^(1-alpha/2)*t(chol(Sigma))%*%matrix(rnorm(n*d,sd=1/sqrt(n)),nrow=d)
  BM <- apply(U,1,cumsum)
  ## Adding the drift, BM and CPP
  X <- drift+BM+sim_CPP(n,avg_no_jumps)
  X
}


###########################
### Plotting some paths ###
###########################

X <- sim_X(1000)
plot(X[,1], type = "l", col = "red", ylim = c(min(X[,c(1,7,11)]),max(X[,c(1,7,11)])))
lines(X[,7], col = "blue")
lines(X[,11], col = "darkgreen")