######################
### Required files ###
######################

## The code below uses functions from tree_learning.R

## Loading data - make sure to set working directory or modify the path
log_prices <- read.csv("data_all.csv")


#############################
### Selecting some stocks ###
#############################

## Selecting certain stocks
selected_tickers <- c("AAPL", "MSFT",
                      "OXY",  "CVX",
                      "UPS",  "FDX",
                      "KO",   "PEP",
                      "HD",   "TGT",
                      "JPM",  "WFC")

selected_data <- log_prices[selected_tickers]


#####################
### Tree learning ###
#####################

## Calculating increments
selected_incs <- diff(data.matrix(selected_data))

## Choosing q
q <- 0.94

## Estimating and plotting the tree
tree <- learn_tree(selected_incs, q)

plot(tree, asp = 0, vertex.label = selected_tickers)


###################
### Subsampling ###
###################

## Some setup
n <- nrow(selected_incs)
d <- ncol(selected_incs)
sub_size <- round(n / 2) #Size of random subsamples
N <- 10000 #Number of subsampling iterations
edge_counts <- matrix(0, nrow = d, ncol = d)
q_range <- c(0.925, 0.95) #q is sampled uniformly in this range

## Function for updating edge_counts
update_counts <- function(incs, q) {
  tree <- learn_tree(incs, q)
  edges <- apply(as_edgelist(tree),1,toString)
  for (i in 1:(d-1)) {
    for (j in (i+1):d) {
      edge_counts[i,j] <<- edge_counts[i,j] + paste(i, j, sep = ", ") %in% edges
    }
  }
}

## Doing the subsampling
## (Takes some time to run)
index_mat <- replicate(N, sample(1:n, size = sub_size, replace = FALSE))
q_vec <- runif(N, min = q_range[1], max = q_range[2])

sapply(1:N, function(j) update_counts(selected_incs[index_mat[,j],], q_vec[j]))

## Plotting fully connected graph with edge widths proportional to edge_counts
weights_fully <- matrix(nrow=0,ncol=3)
for (i in 1:(d-1)) {
  for (j in (i+1):d) {
    weights_fully <- rbind(weights_fully, c(i, j, edge_counts[i,j]))
  }
}
graph_fully <- graph.data.frame(weights_fully[,1:2],directed=F)
E(graph_fully)$weight <- weights_fully[,3]
plot(graph_fully,edge.width = weights_fully[,3] / (N / 5), vertex.label = selected_tickers)

## Calculating and plotting minimum spanning tree based on edge_counts
weights_mst <- matrix(nrow=0,ncol=3)
for (i in 1:(d-1)) {
  for (j in (i+1):d) {
    weights_mst <- rbind(weights_mst, c(i, j, -edge_counts[i,j]))
  }
}
graph_mst <- graph.data.frame(weights_mst[,1:2],directed=F)
E(graph_mst)$weight <- weights_mst[,3]
mst <- minimum.spanning.tree(graph_mst)
plot(mst, asp = 0, vertex.label = selected_tickers)


###########################
### Plotting -log(\chi) ###
###########################

## Some setup
q_range <- seq(from = 0.8, to = 0.99, by = 0.01)
R <- length(q_range)
chi_array <- matrix(nrow = R, ncol = d - 1)

## Estimating chi_{1j} values for j>=2 and for all q in q_range
for (k in 1:R) {
  for (j in 1:(d - 1)) {
    chi_array[k,j] <- chi_est(selected_incs[,1], selected_incs[,j + 1], q_range[k])
  }
}

weights_array <-  -log(chi_array)

## Plotting the weights -log(chi_{1j}) as a function of q
layout(matrix(c(1,2),nrow=1), width=c(4,1)) 
par(mar=c(5,4,4,0))
matplot(q_range, weights_array, type = "l", xlab = "q", ylab = "Weights", lty = 1)
par(mar=c(5,0,4,2))
plot(c(0,1),type="n", axes=F, xlab="", ylab="")
legend("right", sprintf("%s", 2:d),
       col = seq_len(d - 1), fill = seq_len(d - 1), cex=0.8)
