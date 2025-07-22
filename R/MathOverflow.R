#' MathOverflow multilayer network data
#'
#' A directed multilayer network of interactions on the stack exchange website
#' Math Overflow, with 3 layers corresponding to users’ answers to questions(a2q),
#' comments to questions(c2q), and comments to answers(c2a). The network includes 23,654
#' users, i.e., the multilayer network consists of 23,654 nodes.
#'
#' @docType data
#' @name MathOverflow
#' @usage data(MathOverflow)
#' @format The \code{MathOverflow} object is a list of sparse matrices representing
#' the adjacency matrices of the MathOverflow interaction network.
#' All matrices share the same node set with consistent node IDs and dimnames.
#'
#' @references  Paranjape, A., Benson, A. R., and Leskovec, J. (2017).
#' \emph{Motifs in temporal networks.}, \emph{In Proceedings of the tenth ACM
#' International Conference on Web Search and Data Mining.}\cr
#'
#' @examples
#' data(MathOverflow)
#' mlA <- MathOverflow
#' timing <- system.time({
#' cluster <- mrcoclust(mlA, rank.r = 2, rank.c = 3, kr = 2, kc = 3, q = 4, p = 0.9)
#' })
#' print(timing)

"MathOverflow"
