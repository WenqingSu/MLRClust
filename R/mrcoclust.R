#' Randomized spectral co-clustering for multilayer stochastic co-block models
#'
#' Randomized spectral co-clustering for multilayer directed networks. Randomized
#' clustering uses both the random sampling strategy and the random projection
#' strategy. Can deal with very large networks.
#'
#' This function computes the clusters of multilayer directed networks using
#' randomized spectral co-clustering algorithms. Random sampling is first performed
#' on the adjacency matrices, then the random projection-based eigendecomposition
#' is performed on the aggregated matrix. The k-means is then performed on the
#' randomized eigenvectors.
#'
#'
#' @param mlA The adjacency matrices list of a multilayer directed network consists of "\code{dgCMatrix}" types.
#' @param rank.r The target row rank of the low-rank decomposition.
#' @param rank.c The target column rank of the low-rank decomposition.
#' @param kr The number of row clusters.
#' @param kc The number of column clusters.
#' @param q The power parameter. It need to be a positive integer number. Default is 2.
#' @param p The sampling probability. Should be between 0 and 1. Default is 0.7.
#' @param iter.max Maximum number of iterations in the \code{\link[stats]{kmeans}}. Default is 50.
#' @param nstart The number of random sets in \code{\link[stats]{kmeans}}. Default is 10.
#' @param nthread Maximum number of threads for specific computations which could be implemented in parallel. Default is 1.
#'
#' @return \item{row.cluster}{The row cluster vector (from \code{1:kr}) with the numbers indicating which row cluster each node is assigned.}
#'         \item{col.cluster}{The column cluster vector (from \code{1:kc}) with the numbers indicating which column cluster each node is assigned.}
#' @export mrcoclust
#' @examples
#' n <- 600
#' kr <- 2
#' kc <- 2
#' row.label.true <- sample(rep(1:kr, each = n/kr), n)
#' col.label.true <- sample(rep(1:kc, each = n/kc), n)
#' probmat1 <- matrix(0.05, kr, kc) + diag(c(0.1, 0.05))
#' probmat1[1,2] <- 0.08
#' probmat2 <- matrix(0.01, kr, kc) + diag(c(0.05, 0.02))
#' probmat2[2,1] <- 0.1
#' probmats <- list(probmat1, probmat2, probmat2)
#' mlA <- mln_generator(row.label.true, col.label.true, probmats, directed = TRUE)
#' rank.r <- 2
#' rank.c <- 2
#' q <- 2
#' p <- 0.7
#' cc <- mrcoclust(mlA, rank.r, rank.c, kr, kc, q = q, p = p)
#'

mrcoclust <- function(mlA, rank.r, rank.c, kr, kc, q = 2, p = 0.7, iter.max = 50, nstart = 10, nthread = 1) {
  mlA_s <- ml.rsample(mlA, p = p, sym = FALSE)
  n <- nrow(mlA_s[[1]])
  L <- length(mlA_s)

  # Ml.r = AlAl'/np^2 - Dl/np^2, M = mean_l(Ml.r)
  # K.r = [M * Omega, (MM') * M * Omega, ..., (MM')^q * M * Omega]
  K.r <- matrix(numeric(0), nrow = n, ncol = 0)
  # out degree
  mlDout <- lapply(mlA_s, function(x) {
    row_sums <- rowSums(x) / (L * n *  p ** 2)
    Diagonal(x = row_sums)
  })
  Omega <- matrix(rnorm(n * rank.r), nrow = n, ncol = rank.r)
  mlAcoord <- lapply(mlA_s, function(A) sparse_matrix_coords(A, nthread))
  hq <- 2 * q + 1 # Highest power
  for (i in 1:hq) {
    tmp <- lapply(mlAcoord, function(Acoord) krylov_space_part(Acoord, Omega, scale = p * sqrt(n * L), q = 1, nthread = nthread))
    tmp <- mapply("-", tmp, lapply(mlDout, function(D) D %*% Omega), SIMPLIFY = FALSE)
    Omega <- as.matrix(Reduce("+", tmp))
    if (i %% 2 == 1) {
      K.r <- cbind(K.r, Omega)
    }
  }
  # Orthogonalize K using QR decomposition: K = Q * R
  qr_Q_inplace(K.r)
  Q.r <- K.r

  # Obtain the smaller matrix B := Q' * M * Q
  B.r <- lapply(mlAcoord, function(Acoord) krylov_space_part(Acoord, Q.r, scale = p * sqrt(n * L), q = 1, nthread = nthread))
  B.r <- mapply("-", B.r, lapply(mlDout, function(D) D %*% Q.r), SIMPLIFY = FALSE)
  B.r <- crossprod(Q.r, as.matrix(Reduce("+", B.r)))

  # Compute the eigenvalue decomposition of B
  fit.r <- eigen(B.r, symmetric = TRUE)
  o <- order(abs(fit.r$values), decreasing = T)
  vectors.r <- Q.r %*% fit.r$vectors[, o[0:rank.r]]

  # Ml.c = Al'Al/np^2 - Dl/np^2, M = mean_l(Ml.c)
  # K.c = [M * Omega, (MM') * M * Omega, ..., (MM')^q * M * Omega]
  K.c <- matrix(numeric(0), nrow = n, ncol = 0)
  # in degree
  mlDin <- lapply(mlA_s, function(x) {
    col_sums <- colSums(x) / (L * n *  p ** 2)
    Diagonal(x = col_sums)
  })
  Omega <- matrix(rnorm(n * rank.c), nrow = n, ncol = rank.c)
  mlA_s <- lapply(mlA_s, t)
  mlAcoord <- lapply(mlA_s, function(A) sparse_matrix_coords(A, nthread))
  hq <- 2 * q + 1 # Highest power
  for (i in 1:hq) {
    tmp <- lapply(mlAcoord, function(Acoord) krylov_space_part(Acoord, Omega, scale = p * sqrt(n * L), q = 1, nthread = nthread))
    tmp <- mapply("-", tmp, lapply(mlDin, function(D) D %*% Omega), SIMPLIFY = FALSE)
    Omega <- as.matrix(Reduce("+", tmp))
    if (i %% 2 == 1) {
      K.c <- cbind(K.c, Omega)
    }
  }
  # Orthogonalize K using QR decomposition: K = Q * R
  qr_Q_inplace(K.c)
  Q.c <- K.c

  # Obtain the smaller matrix B := Q' * M * Q
  B.c <- lapply(mlAcoord, function(Acoord) krylov_space_part(Acoord, Q.c, scale = p * sqrt(n * L), q = 1, nthread = nthread))
  B.c <- mapply("-", B.c, lapply(mlDin, function(D) D %*% Q.c), SIMPLIFY = FALSE)
  B.c <- crossprod(Q.c, as.matrix(Reduce("+", B.c)))

  # Compute the eigenvalue decomposition of B
  fit.c <- eigen(B.c, symmetric = TRUE)
  o <- order(abs(fit.c$values), decreasing=T)
  vectors.c <- Q.c %*% fit.c$vectors[, o[0:rank.c]]

  # clustering
  fit.r <- kmeans(vectors.r, kr, iter.max = 50, nstart = 10)
  cluster.r <- fit.r$cluster
  fit.c <- kmeans(vectors.c, kc, iter.max = 50, nstart = 10)
  cluster.c <- fit.c$cluster

  list(row.cluster = cluster.r, col.cluster = cluster.c)
}
