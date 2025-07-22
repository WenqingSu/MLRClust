#' Generate multilayer networks from multilayer stochastic block models
#' and multilayer stochastic co-block models.
#'
#' Generate adjacency matrices of multilayer undirected networks from
#' multilayer stochastic block models and multilayer directed networks from
#' multilayer stochastic co-block models.
#'
#' This function generates multilayer undirected (directed) networks using
#' multilayer stochastic block (co-block) models.In the multilayer stochastic
#' block (co-block) models, all layers share common row and column communities
#' but with possibly different edge densities. The layer-wise networks are
#' generated independently from the stochastic block (co-block) model.


#' @param directed If \code{TRUE}, the multilayer directed networks is generated from the multilayer stochastic co-block model.
#'                 If \code{FALSE}, the multilayer undirected networks is generated from the multilayer stochastic block model. Default is \code{FALSE}.
#' @param row.label The row community vector (from \code{1:kr}) with the numbers indicating which row community each node is assigned.
#' @param col.label The column community vector (from \code{1:kc}) with the numbers indicating which column community each node is assigned.
#'                  \code{row.label} and \code{col.label} must have the same length but not necessarily indicate the same communities and community numbers.
#'                  Only required when \code{directed = TRUE}.
#' @param probmat List of link probability matrices with dimension \code{c(kr, kc)}.
#'                \code{probmat[l][i,j]} is proportional to the probability of an edge from nodes in row community \code{i} to nodes in column community \code{j} in layer \code{l}.
#' @param sparse Logical indicating whether to return sparse matrices (\code{dgCMatrix}) for network adjacency matrices.
#'               Default is \code{TRUE}.
#' @return A list of adjacency matrices representing the generated multilayer networks.
#' @references W. Su, X. Guo, X. Chang and Y. Yang. (2024)
#' \emph{Spectral co-clustering in multi-layer directed networks},
#' \emph{Computational Statistics & Data Analysis, Vol. 198, 107987}\cr
#' \url{https://doi.org/10.1016/j.csda.2024.107987}\cr
#'
#' @export mln_generator
#' @examples
#' # The multilayer stochastic block model
#' n <- 300
#' kr <- 2
#' row.label <- sample(rep(1:kr, each = n/kr), n)
#' probmat1 <- matrix(0.05, kr, kr) + diag(c(0.1, 0.05))
#' probmat2 <- matrix(0.01, kr, kr) + diag(c(0.05, 0.02))
#' probmats <- list(probmat1, probmat2, probmat2)
#' mln1 <- mln_generator(row.label = row.label, probmat = probmats, directed = FALSE)
#'
#' # The multilayer stochastic co-block model
#' n <- 300
#' kr <- 2
#' kc <- 3
#' row.label <- sample(rep(1:kr, each = n/kr), n)
#' col.label <- sample(rep(1:kc, each = n/kc), n)
#' probmat1 <- matrix(0.05, kr, kc)
#' probmat1[1,1] <- 0.1
#' probmat1[2,2] <- 0.1
#' probmat2 <- matrix(0.01, kr, kc)
#' probmat1[1,1] <- 0.05
#' probmat1[2,3] <- 0.03
#' probmats <- list(probmat1, probmat2, probmat2)
#' mln2 <- mln_generator(row.label, col.label, probmats, directed = TRUE)
#'

mln_generator <- function(row.label, col.label, probmat, directed = FALSE, sparse = TRUE) {
  n_nodes <- length(row.label)
  layers <- length(probmat)
  adj <- matrix(0, nrow = n_nodes, ncol = n_nodes)
  adjs <- vector("list", length = layers)
  if (directed) {
    if (missing(col.label)) {
      stop("The 'col.label' argument is required but not provided.")
    }
    for (l in 1:layers) {
      p_matrix <- probmat[[l]][row.label, , drop = FALSE]
      p_matrix <- p_matrix[, col.label, drop = FALSE]
      adj <- array(rbinom(n = prod(dim(p_matrix)), size = 1, prob = c(p_matrix)), dim = dim(p_matrix))
      diag(adj) <- 0
      if (sparse)
        adjs[[l]] <- as(Matrix(adj, sparse = TRUE), "dgCMatrix")
      else
        adjs[[l]] <- adj
    }
  }
  else {
    for (l in 1:layers) {
      p_matrix <- probmat[[l]][row.label, , drop = FALSE]
      p_matrix <- p_matrix[, row.label, drop = FALSE]
      # Generate upper triangular matrix
      upper_tri_indices <- which(upper.tri(p_matrix), arr.ind = TRUE)
      adj[upper_tri_indices] <- rbinom(n = length(upper_tri_indices)/2, size = 1, prob = p_matrix[upper_tri_indices])
      # Mirror upper to lower triangular part
      adj[lower.tri(adj)] <- t(adj)[lower.tri(adj)]
      diag(adj) <- 0
      if (sparse)
        adjs[[l]] <- as(Matrix(adj, sparse = TRUE), "dgCMatrix")
      else
        adjs[[l]] <- adj
    }
  }
  return(adjs)
}
