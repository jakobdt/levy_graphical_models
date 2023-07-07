######################
### Required files ###
######################

## The code below uses functions from simulation.R and tree_learning.R


########################################################
### Monte Carlo estimation of recovery probabilities ###
########################################################

## Function for estimating probability of recovering the tree.
## K is the number of simulations to base each probability estimate on,
## qvec is the values of q to use for the structure learning, and n is
## the number of times the process is sampled (i.e. sampling frequency is 1/n)
est_probs <- function(K, qvec, n) {
  X <- replicate(K, sim_X(n))
  d <- dim(X)[2]
  results <- vector(length = length(qvec))
  for (i in 1:length(qvec)) {
    matches <- sapply(1:K, function(k) length(E(graph.intersection(Tree, learn_tree(X[,,k], qvec[i])))) == d - 1)
    results[i] <- mean(matches)
  }
  results
}

## Setting the values to be used
K <- 1000
q_range <- seq(.5, .99, .01)
n_range <- c(200, 500, 1000, 2000, 5000)

## Estimation of probabilities
results <- matrix(nrow = length(q_range), ncol = length(n_range))
for (j in 1:length(n_range)) {
  results[,j] <- est_probs(K, q_range, n_range[j])
}

## Plotting probabilities as a function of q
N <- length(n_range)
layout(matrix(c(1,2),nrow=1), width=c(4,1)) 
par(mar=c(5,4,4,0))
matplot(q_range, results, ylim = c(0,1), type = "l", xlab = "q", ylab = "Recovery probability")
par(mar=c(5,0,4,2))
plot(c(0,1),type="n", axes=F, xlab="", ylab="")
legend("right", c("200", "500", "1000", "2000", "5000"),
       col = seq_len(N), fill = seq_len(N), cex=0.8)