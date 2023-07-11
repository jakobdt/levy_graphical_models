################
### Packages ###
################

require(igraph)


################################
### Functions for estimation ###
################################

## Function for computing increments of a process X
increments <- function(X) {
  n <- length(X)
  X-c(0,X[1:(n-1)])
}

## Applies empirical CDF of increments to the increments
ecdf_transform <- function(incs) {
  n <- length(incs)
  rank(incs) / (n + 1)
}

## Function for estimating chi based on positive quadrant
positive_quadrant <- function(incs1, incs2, q) {
  n <- length(incs1)
  Z1 <- ecdf_transform(incs1)
  Z2 <- ecdf_transform(incs2)
  sum((Z1 > q) * (Z2 > q)) / (n * (1 - q))
}

## Function for estimation chi based on all quadrants
chi_est <- function(incs1, incs2, q) {
  res <- 0
  for (s1 in c(-1,1)) {
    for (s2 in c(-1,1)) {
      res <- res + positive_quadrant(s1 * incs1, s2 * incs2, q)
    }
  }
  res
}

## Function for learning the tree
learn_tree <- function(incs, q) {
  d = ncol(incs)
  chi <- matrix(nrow = d, ncol = d)
  for (i in 1:(d - 1)) {
    for (j in (i + 1):d) {
      chi[i,j] <- chi_est(incs[,i], incs[,j], q)
    }
  }
  weight_matrix <- matrix(nrow=0,ncol=3)
  for (i in 1:(d-1)) {
    for (j in (i+1):d) {
      weight_matrix <- rbind(weight_matrix,c(i,j,-log(chi[i,j])))
    }
  }
  weighted_graph <- graph.data.frame(weight_matrix[,1:2],directed=F)
  E(weighted_graph)$weight <- weight_matrix[,3]
  mst <- minimum.spanning.tree(weighted_graph)
  mst
}


######################################################
### Example (requires functions from simulation.R) ###
######################################################

X <- sim_X(5000)
learn_tree(increments(X), 0.9)