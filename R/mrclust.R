#' Randomized spectral clustering for multilayer stochastic block models
#'
#' Randomized spectral clustering for multilayer undirected networks. Randomized 
#' clustering uses both the random sampling strategy and the random projection 
#' strategy. Can deal with very large networks.
#'
#' This function computes the clusters of multilayer undirected networks using
#' randomized spectral clustering algorithms. Random sampling is first performed 
#' on the adjacency matrices, then the random projection-based eigendecomposition 
#' is performed on the aggregated matrix. The k-means is then performed on the 
#' randomized eigenvectors.
#'
#'
#' @param mlA The adjacency matrices list of a multilayer undirected network consists of "\code{dgCMatrix}" types.
#' @param rank The target rank of the low-rank decomposition.
#' @param k The number of clusters.
#' @param q The power parameter. It need to be a positive integer number. Default is 2.
#' @param p The sampling probability. Should be between 0 and 1. Default is 0.7.
#' @param iter.max Maximum number of iterations in the \code{\link[stats]{kmeans}}. Default is 50.
#' @param nstart The number of random sets in \code{\link[stats]{kmeans}}. Default is 10.
#' @param nthread Maximum number of threads for specific computations which could be implemented in parallel. Default is 1.
#'
#' @return \item{cluster}{The cluster vector (from \code{1:k}) with the numbers indicating which
#'              cluster each node is assigned.}
#' @export mrclust
#' @examples
#' n <- 500
#' k <- 2
#' label.true <- sample(rep(1:k, each = n/k), n)
#' probmat1 <- matrix(0.05, k, k) + diag(c(0.1, 0.05))
#' probmat2 <- matrix(0.01, k, k) + diag(c(0.05, 0.02))
#' probmats <- list(probmat1, probmat2, probmat2)
#' mlA <- mln_generator(directed = FALSE, label.true, probmat = probmats)
#' rank <- 2
#' q <- 2
#' p <- 0.7
#' mrclust(mlA, rank, k = k, q = q, p = p)
#'

mrclust <- function(mlA, rank, k, q = 2, p = 0.7, iter.max = 50, nstart = 10, nthread = 1) {
  mlA_s <- ml.rsample(mlA, p = p, sym = TRUE)
  n <- nrow(mlA_s[[1]])
  L <- length(mlA_s)
  
  # Ml = AlAl'/np^2 - Dl/np^2, M = mean_l(Ml)
  # K = [M * Omega, (MM') * M * Omega, ..., (MM')^q * M * Omega]
  K <- matrix(numeric(0), nrow = n, ncol = 0)
  # out degree
  mlDout <- lapply(mlA_s, function(x) {
    row_sums <- rowSums(x) / (L * n *  p ** 2)
    Diagonal(x = row_sums)
  }) 
  Omega <- matrix(rnorm(n * rank), nrow = n, ncol = rank)
  mlAcoord <- lapply(mlA_s, function(A) sparse_matrix_coords(A, nthread))
  hq <- 2 * q + 1 # Highest power
  for (i in 1:hq) {
    tmp <- lapply(mlAcoord, function(Acoord) krylov_space_part(Acoord, Omega, scale = p * sqrt(n * L), q = 1, nthread = nthread))
    tmp <- mapply("-", tmp, lapply(mlDout, function(D) D %*% Omega), SIMPLIFY = FALSE)
    Omega <- as.matrix(Reduce("+", tmp))
    if (i %% 2 == 1) {
      K <- cbind(K, Omega)
    }
  }
  # Orthogonalize K using QR decomposition: K = Q * R
  qr_Q_inplace(K)
  Q <- K
  
  # Obtain the smaller matrix B := Q' * M * Q
  B <- lapply(mlAcoord, function(Acoord) krylov_space_part(Acoord, Q, scale = p * sqrt(n * L), q = 1, nthread = nthread))
  B <- mapply("-", B, lapply(mlDout, function(D) D %*% Q), SIMPLIFY = FALSE)
  B <- crossprod(Q, as.matrix(Reduce("+", B)))
  
  # Compute the eigenvalue decomposition of B
  fit <- eigen(B, symmetric = TRUE)
  o <- order(abs(fit$values), decreasing = T)
  vectors <- Q %*% fit$vectors[, o[0:rank]]
  
  fit <- kmeans(vectors, k, iter.max = 50, nstart = 10)
  cluster <- fit$cluster
  
  cluster
}