ml.rsample <- function(mlA, p, sym = TRUE) {
  L <- length(mlA)
  if (sym) {
    for (l in 1:L) {
      mlA[[l]] <- rsample_sym(mlA[[l]], p) # strictly lower triangular matrix
      mlA[[l]] <- mlA[[l]] + t(mlA[[l]])}
  }
  else {
    for (l in 1:L)
      mlA[[l]] <- rsample(mlA[[l]], p)
  }
  mlA
}
