min_mis_error <- function(true_label, kmeans_labels, k_clusters) {
  location <- vector("list", k_clusters)
  for (i in 1:k_clusters) {
    location[[i]] <- which(kmeans_labels == i)
  }

  error_min <- 1
  best_p <- NULL

  permutations_ <- permutations(k_clusters, k_clusters)

  for (cp in 1:nrow(permutations_)) {
    perm <- permutations_[cp, ]
    for (j in 1:k_clusters) {
      kmeans_labels[location[[j]]] <- perm[j]
    }
    error <- sum(kmeans_labels != true_label) / length(true_label)
    if (error <= error_min) {
      error_min <- error
      best_p <- perm
    }
  }
  return(error_min)
}


dsos <- function(mlA, rank.t, k){
  n <- nrow(mlA[[1]])
  L <- length(mlA)

  mlA_square <- lapply(mlA, function(A) {
    Asquare <- tcrossprod(A) / (n * L)
  })

  M <- Reduce("+", mlA_square)
  diag(M) <- 0
  M

  vectors <- partial_eigen(M, n = rank.t)

  fit <- kmeans(vectors$vectors, k, iter.max = 50, nstart = 10)
  cluster <- fit$cluster

  return(cluster)
}
