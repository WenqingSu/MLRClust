% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrclust.R
\name{mrclust}
\alias{mrclust}
\title{Randomized spectral clustering for multilayer stochastic block models}
\usage{
mrclust(mlA, rank, k, q = 2, p = 0.7, iter.max = 50, nstart = 10, nthread = 1)
}
\arguments{
\item{mlA}{The adjacency matrices list of a multilayer undirected network consists of "\code{dgCMatrix}" types.}

\item{rank}{The target rank of the low-rank decomposition.}

\item{k}{The number of clusters.}

\item{q}{The power parameter. It need to be a positive integer number. Default is 2.}

\item{p}{The sampling probability. Should be between 0 and 1. Default is 0.7.}

\item{iter.max}{Maximum number of iterations in the \code{\link[stats]{kmeans}}. Default is 50.}

\item{nstart}{The number of random sets in \code{\link[stats]{kmeans}}. Default is 10.}

\item{nthread}{Maximum number of threads for specific computations which could be implemented in parallel. Default is 1.}
}
\value{
\item{cluster}{The cluster vector (from \code{1:k}) with the numbers indicating which
             cluster each node is assigned.}
}
\description{
Randomized spectral clustering for multilayer undirected networks. Randomized 
clustering uses both the random sampling strategy and the random projection 
strategy. Can deal with very large networks.
}
\details{
This function computes the clusters of multilayer undirected networks using
randomized spectral clustering algorithms. Random sampling is first performed 
on the adjacency matrices, then the random projection-based eigendecomposition 
is performed on the aggregated matrix. The k-means is then performed on the 
randomized eigenvectors.
}
\examples{
n <- 500
k <- 2
label.true <- sample(rep(1:k, each = n/k), n)
probmat1 <- matrix(0.05, k, k) + diag(c(0.1, 0.05))
probmat2 <- matrix(0.01, k, k) + diag(c(0.05, 0.02))
probmats <- list(probmat1, probmat2, probmat2)
mlA <- mln_generator(directed = FALSE, label.true, probmat = probmats)
rank <- 2
q <- 2
p <- 0.7
mrclust(mlA, rank, k = k, q = q, p = p)

}
