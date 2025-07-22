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
#' @param rank.t The target rank of the low-rank decomposition.
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
#' # example 1
#' n <- 500
#' k <- 2
#' label.true <- sample(rep(1:k, each = n/k), n)
#' probmat1 <- matrix(0.05, k, k) + diag(c(0.1, 0.05))
#' probmat2 <- matrix(0.01, k, k) + diag(c(0.05, 0.02))
#' probmats <- list(probmat1, probmat2, probmat2)
#' mlA <- mln_generator(row.label = label.true, probmat = probmats, directed = FALSE)
#' rank.t <- 2
#' q <- 2
#' p <- 0.7
#' mrclust(mlA, rank.t, k = k, q = q, p = p)
#'
#'
#' # example 2
#' # The effect of power parameter.
#' U <- matrix(c(1/2, 1/2, -sqrt(2)/2, 1/2, 1/2, sqrt(2)/2, sqrt(2)/2, -sqrt(2)/2, 0), nrow=3, byrow=TRUE)
#' Lambda1 <- diag(c(1.5, 0.2, 0.4))
#' Lambda2 <- diag(c(1.5, 0.2, -0.4))
#' B_1 <- U %*% Lambda1 %*% t(U)
#' B_2 <- U %*% Lambda2 %*% t(U)
#' layer <- 20
#' density <- 0.1
#' probmats <- lapply(1:layer, function(l) {if (l <= floor(layer / 2)) B_1 * density else B_2 * density})
#' k <- 3
#' p <- 0.7 # sampling probability
#' nodes_list <- seq(400, 2000, 200)
#' mr.result <- matrix(0, length(nodes_list), 4) # misclassification rates
#' for (i in 1:length(nodes_list)) {
#'   n_nodes <- nodes_list[i]
#'   size_communities <- c(0.3, 0.4, 0.3) * n_nodes
#'   true_label <- sample(rep(seq_len(k), times=size_communities), n_nodes)
#'     for (j in 1:5){
#'       mlA <- mln_generator(row.label = true_label, probmat = probmats, directed = FALSE)
#'       # randomized spectral clustering q = 2
#'       cluster.r.q2 <- mrclust(mlA, rank.t = k, k = k, q = 2, p = p)
#'       mr.result[i, 1] <- mr.result[i, 1] + min_mis_error(true_label, cluster.r.q2, k_clusters = k) / 5
#'       # randomized spectral clustering q = 4
#'       cluster.r.q4 <- mrclust(mlA, rank.t = k, k = k, q = 4, p = p)
#'       mr.result[i, 2] <- mr.result[i, 2] + min_mis_error(true_label, cluster.r.q4, k_clusters = k) / 5
#'       # randomized spectral clustering q = 6
#'       cluster.r.q6 <- mrclust(mlA, rank.t = k, k = k, q = 6, p = p)
#'       mr.result[i, 3] <- mr.result[i, 3] + min_mis_error(true_label, cluster.r.q6, k_clusters = k) / 5
#'       # original spectral clustering
#'       cluster.sc <- dsos(mlA, rank.t = k, k = k)
#'       mr.result[i, 4] <- mr.result[i, 4] + min_mis_error(true_label, cluster.sc, k_clusters = k) / 5
#'   }
#' }
#' line_labels <- c("RSC q = 2", "RSC q = 4", "RSC q = 6", "SC")
#' matplot(nodes_list, mr.result, type = "l", lty = 1, col = 1:4, xlab = "Number of nodes (n)", ylab = "Misclassification rate")
#' legend("topright", legend = line_labels, col = 1:4, lty = 1)

mrclust <- function(mlA, rank.t, k, q = 2, p = 0.7, iter.max = 50, nstart = 10, nthread = 1) {
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
  Omega <- matrix(rnorm(n * rank.t), nrow = n, ncol = rank.t)
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
  vectors <- Q %*% fit$vectors[, o[0:rank.t]]

  fit <- kmeans(vectors, k, iter.max = 50, nstart = 10)
  cluster <- fit$cluster

  cluster
}
